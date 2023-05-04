import os
import re
from glob import glob
import subprocess as sp
import tempfile

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import click


BED_FILE_ENDING = '.bed'
FASTA_FILE_ENDING = '.fa'


def _parse_fasta(buffer):
    seqs = {}
    seq_id = next(buffer).strip().split(' ')[0][1:]
    seq = []
    for line in buffer:
        if not line.startswith('>'):
            seq.append(line.strip())
        else:
            seqs[seq_id] = ''.join(seq)
            seq_id = line.strip().split(' ')[0][1:]
            seq = []
    else:
        seqs[seq_id] = ''.join(seq)
    return seqs


def read_fasta(fasta_fn):
    with open(fasta_fn) as f:
        seqs = _parse_fasta(f)
    return seqs


def read_orthofinder_n0(n0_fn, min_og_size=2):
    og = pd.read_csv(n0_fn, sep='\t', index_col='HOG')
    og = og.iloc[:, 2:]
    og = og[(~og.isnull()).sum(1) >= min_og_size]
    return og


def read_intron_bed(intron_bed_fn):
    intron_data = pd.read_csv(
        intron_bed_fn, sep='[\t|]',
        engine='python',
        names=[
            'chrom', 'start', 'end', 'prot_id', 
            'prot_position', 'prot_frame', 'score', 'strand'
        ],
        usecols=[
            'chrom', 'start', 'end', 'prot_id',
            'prot_position', 'prot_frame', 'strand'
        ],
        dtype={
            'chrom': str, 'start': int, 'end': int,
            'prot_id': str, 'prot_position': int,
            'prot_frame': int, 'strand': str
        }
    )
    return intron_data


def get_species_intron_data(intron_bed_fns):
    species_intron_data = {}
    for species_name, bed_fn in intron_bed_fns.items():
        species_intron_data[species_name] = read_intron_bed(bed_fn)
    return species_intron_data


def get_species_protein_data(protein_fasta_fns):
    species_protein_data = {}
    for species_name, fasta_fn in protein_fasta_fns.items():
        species_protein_data[species_name] = read_fasta(fasta_fn)
    return species_protein_data


def find_ogs_with_introns(og, intron_data, prot_data,
                          min_species_with_introns=2,
                          skip_very_large_og=True):
    if skip_very_large_og:
        n_species = og.shape[1]
        skip_size = n_species * 2
    else:
        skip_size = np.inf
    for og_id, group in og.iterrows():
        group = group.dropna().str.split(', ').explode()
        og_introns = []
        og_intron_species = set()
        for species, prot_id in group.items():
            species_intron_data = intron_data[species]
            i = species_intron_data.query('prot_id == @prot_id')
            if len(i):
                i = i.assign(species_name=species, og_id=og_id)
                og_introns.append(i)
                og_intron_species.add(species)
        else:
            if len(og_intron_species) >= min_species_with_introns:
                og_introns = pd.concat(og_introns)
                og_seqs = {}
                og_prot_ids = og_introns[['species_name', 'prot_id']].drop_duplicates()
                if len(og_prot_ids) < skip_size:
                    for _, species, prot_id in og_prot_ids.itertuples():
                        og_seqs[prot_id] = prot_data[species][prot_id]
                    yield og_introns, og_seqs


def sanitise_protein_seqs(seqs):
    cleaned = {
        seq_id: re.sub('[^ACDEFGHIKLMNPQRSTVWY]', 'X', seq)
        for seq_id, seq in seqs.items()
    }
    return cleaned


def run_mafft(seqs, max_iter=1000, fast_method_threshold=200):
    tmp_fh, tmp_fn = tempfile.mkstemp()
    if len(seqs) <= fast_method_threshold:
        cmd = "mafft --localpair --maxiterate {max_iter} {tmp_fn}"
    else:
        cmd = "mafft --retree 2 --maxiterate {max_iter} {tmp_fn}"
    cmd = cmd.format(max_iter=max_iter, tmp_fn=tmp_fn)
    with open(tmp_fn, 'w') as tmp:
        for seq_id, seq in sanitise_protein_seqs(seqs).items():
            tmp.write(f">{seq_id}\n{seq}\n")
    with sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE) as proc:
        in_data = []
        dout, derr = proc.communicate()
    if proc.returncode:
        print(derr.decode())
        raise IOError()
    dout = dout.decode()
    os.close(tmp_fh)
    os.remove(tmp_fn)
    return _parse_fasta(iter(dout.split('\n')))


def aligned_seq_to_aligned_pos(aligned_seq):
    pos = np.array([s != '-' for s in aligned_seq]).astype(int)
    pos_cs = np.cumsum(pos)
    pos_cs = np.ma.masked_array(pos_cs, mask=pos == 0)
    return pos_cs


def get_msa_column_label(aligned_pos, prot_pos, prot_frame):
    column = np.searchsorted(aligned_pos, prot_pos)
    labels = []
    for col, frame in zip(column, prot_frame):
        labels.append(f'col{col}.frame{frame}')
    return labels


def _groupby_transform_msa_col_label(group, *, all_seq_pos):
    seq_pos = all_seq_pos[group.name]
    col_ids = get_msa_column_label(seq_pos, group.prot_position, group.prot_frame)
    return group.assign(col_id = col_ids)


def label_orthogroup_introns_by_msa_column(introns, prot_seqs, min_species_per_col_id=2):
    aligned_seqs = run_mafft(prot_seqs)
    if not aligned_seqs:
        # mafft failed for some reason
        # will fix later
        return pd.DataFrame()
    aligned_seq_pos = {seq_id: aligned_seq_to_aligned_pos(seq)
                       for seq_id, seq in aligned_seqs.items()}
    introns_labelled = introns.groupby('prot_id', as_index=False).apply(
        _groupby_transform_msa_col_label, all_seq_pos=aligned_seq_pos
    )
    introns_labelled = introns_labelled.reset_index(drop=True)
    introns_labelled = introns_labelled.groupby('col_id').filter(
        lambda group: len(group.drop_duplicates('species_name')) >= min_species_per_col_id
    )
    return introns_labelled


def intron_conservation_analysis(og_fn, intron_bed_fns, protein_fasta_fns,
                                 skip_very_large_og=True, n_jobs=1):
    og = read_orthofinder_n0(og_fn)
    intron_data = get_species_intron_data(intron_bed_fns)
    prot_data = get_species_protein_data(protein_fasta_fns)
    og_intron_iterator = find_ogs_with_introns(
        og,
        intron_data,
        prot_data,
        skip_very_large_og=skip_very_large_og
    )
    with Parallel(n_jobs=n_jobs) as pool:
        results = pool(
            delayed(label_orthogroup_introns_by_msa_column)(
                og_introns, og_seqs
            ) for og_introns, og_seqs in og_intron_iterator
        )
    results = [df for df in results if len(df)]
    return pd.concat(results)


def glob_file_type(directory, file_type):
    fn_pattern = os.path.join(
        directory, '*' + file_type
    )
    fns = glob(fn_pattern)
    fns = {
        os.path.split(fn)[1][:-len(file_type)]: fn
        for fn in fns
    }
    return fns


@click.command()
@click.option('-o', '--orthogroup-fn', required=True)
@click.option('-b', '--intron-bed-dir', required=True)
@click.option('-f', '--protein-fasta-dir', required=True)
@click.option('-d', '--output-intron-bed-dir', required=True)
@click.option('--skip-very-large-og/--no-skip-very-large-og', required=False, default=True)
@click.option('-p', '--processes', required=False, default=1)
def main(orthogroup_fn, intron_bed_dir, protein_fasta_dir,
         output_intron_bed_dir, skip_very_large_og, processes):
    intron_bed_fns = glob_file_type(intron_bed_dir, BED_FILE_ENDING)
    protein_fasta_fns = glob_file_type(protein_fasta_dir, FASTA_FILE_ENDING)
    species = list(set(intron_bed_fns).intersection(protein_fasta_fns))
    if not len(species):
        raise IOError('could not find any matching bed and fasta files')
    else:
        intron_bed_fns = {sn: fn for sn, fn in intron_bed_fns.items() if sn in species}
        protein_fasta_fns = {sn: fn for sn, fn in protein_fasta_fns.items() if sn in species}

    if not os.path.exists(output_intron_bed_dir):
        os.mkdir(output_intron_bed_dir)

    results = intron_conservation_analysis(
        orthogroup_fn, intron_bed_fns, protein_fasta_fns,
        skip_very_large_og=skip_very_large_og,
        n_jobs=processes
    )
    # this guarantees that empty categories (species with no conserved introns)
    # are still reported
    results['species_name'] = pd.Categorical(
        results.species_name, categories=species
    )
    # output to one bed file per species
    for species_name, species_introns in results.groupby('species_name'):
        species_output_fn = os.path.join(
            output_intron_bed_dir, f'{species_name}.conserved_introns.bed'
        )
        species_introns = species_introns.sort_values(['chrom', 'start'])
        with open(species_output_fn, 'w') as o:
            for _, i in species_introns.iterrows():
                record = (
                    f'{i.chrom}\t{i.start}\t{i.end}\t'
                    f'{i.prot_id}|{i.prot_position}|{i.prot_frame}|'
                    f'{i.og_id}|{i.col_id}\t.\t{i.strand}\n'
                )
                o.write(record)


if __name__ == '__main__':
    main()