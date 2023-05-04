import re
import os
import warnings
from glob import glob
import tempfile
import subprocess as sp
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import click

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


def read_orthofinder_n0(n0_fn):
    og = pd.read_csv(n0_fn, sep='\t', index_col='HOG').fillna('')
    og = og.iloc[:, 2:]
    return og


def get_species_protein_data(protein_fasta_fns):
    species_protein_data = {}
    for species_name, fasta_fn in protein_fasta_fns.items():
        species_protein_data[species_name] = read_fasta(fasta_fn)
    return species_protein_data


def identify_orthogroups(key_prot_ids, og):
    og_ids = []
    for _, gene_symbol, prot_id, species_name in key_prot_ids.itertuples():
        try:
            og_id = og.query(f'{species_name}.str.contains("{prot_id}")').index[0]
        except IndexError:
            warnings.warn(f'Could not find orthogroup for "{gene_symbol}"')
            og_id = np.nan
        og_ids.append(og_id)
    key_prot_ids['og_id'] = og_ids
    for og_id, group in key_prot_ids.groupby('og_id'):
        yield '_'.join(group.gene_symbol), og_id


def get_og_seqs(og_id, og, seqs):
    og_seqs = {}
    for species, prot_ids in og.loc[og_id].items():
        if len(prot_ids):
            for p_id in prot_ids.split(', '):
                og_seqs[p_id] = seqs[species][p_id]
    return og_seqs

def sanitise_protein_seqs(seqs):
    cleaned = {
        seq_id: re.sub('[^ACDEFGHIKLMNPQRSTVWY]', 'X', seq)
        for seq_id, seq in seqs.items()
    }
    return cleaned


def _parse_hmmscan_tab(lines):
    res = []
    for line in lines.split('\n'):
        line = line.strip()
        if len(line) == 0:
            continue
        if line.startswith('#'):
            continue
        _, _, prot_id, _, _, score, *_ = re.split(r'\s+', line)
        res.append(float(score))
    return res


def build_and_score_hmm(hmm_name, seqs, output_dir, score_perc=2.5,
                        mafft_max_iter=1000, mafft_fast_method_threshold=200):
    tmp_fh, tmp_fn = tempfile.mkstemp()
    output_fn = os.path.join(output_dir, f'{hmm_name}.hmm')
    if len(seqs) <= mafft_fast_method_threshold:
        align_cmd = "mafft --localpair --maxiterate {max_iter} {tmp_fn}"
    else:
        align_cmd = "mafft --retree 2 --maxiterate {max_iter} {tmp_fn}"
    build_cmd = align_cmd + ' | hmmbuild -n {hmm_name} --informat AFA {output_fn} -'
    build_cmd = build_cmd.format(
        max_iter=mafft_max_iter,
        tmp_fn=tmp_fn,
        hmm_name=hmm_name,
        output_fn=output_fn
    )
    with open(tmp_fn, 'w') as tmp:
        for seq_id, seq in sanitise_protein_seqs(seqs).items():
            tmp.write(f">{seq_id}\n{seq}\n")
    build_proc = sp.run(build_cmd, shell=True, check=True,
                        stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    press_cmd = 'hmmpress -f {hmm_fn}'.format(hmm_fn=output_fn)
    press_proc = sp.run(press_cmd, shell=True, check=True,
                        stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    # now run hmm against input seqs to get evalue distribution
    scan_cmd = 'hmmsearch -o /dev/null --tblout /dev/stdout {hmm_fn} {tmp_fn}'
    scan_cmd = scan_cmd.format(hmm_fn=output_fn, tmp_fn=tmp_fn)
    scan_proc = sp.run(scan_cmd, shell=True, check=True,
                       stdout=sp.PIPE, stderr=sp.DEVNULL)
    os.close(tmp_fh)
    os.remove(tmp_fn)
    scores = _parse_hmmscan_tab(scan_proc.stdout.decode())
    return [hmm_name, np.percentile(scores, score_perc)]


def _prepare_hmm_input(n0_fn, key_prot_fn, prot_fasta_fns):
    og = read_orthofinder_n0(n0_fn)
    prot_data = get_species_protein_data(prot_fasta_fns)
    key_prot_table = pd.read_csv(
        key_prot_fn, sep='\t', names=['gene_symbol', 'prot_id', 'species']
    )
    for gene_symbol, og_id in identify_orthogroups(key_prot_table, og):
        og_seqs = get_og_seqs(og_id, og, prot_data)
        yield gene_symbol, og_seqs


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
@click.option('-k', '--key-prot-fn', required=True)
@click.option('-f', '--protein-fasta-dir', required=True)
@click.option('-d', '--output-hmm-dir', required=True)
@click.option('-s', '--score-percentile', default=2.5)
@click.option('-p', '--processes', required=False, default=1)
def main(orthogroup_fn, key_prot_fn, protein_fasta_dir,
         output_hmm_dir, score_percentile, processes):
    protein_fasta_fns = glob_file_type(protein_fasta_dir, '.fa')
    if not os.path.exists(output_hmm_dir):
        os.mkdir(output_hmm_dir)
    og_seq_iterator = _prepare_hmm_input(
        orthogroup_fn, key_prot_fn, protein_fasta_fns
    )
    with Parallel(n_jobs=processes) as pool:
        hmm_info = pool(
            delayed(build_and_score_hmm)(
                hmm_name, seqs, output_hmm_dir, score_percentile,
            ) for hmm_name, seqs in og_seq_iterator
        )
    hmm_info = pd.DataFrame(hmm_info, columns=['hmm_name', 'score_threshold'])
    hmm_info.to_csv(
        os.path.join(output_hmm_dir, 'hmm_scoring_thresholds.tsv'),
        sep='\t',
        float_format='%.3e',
    )


if __name__ == '__main__':
    main()