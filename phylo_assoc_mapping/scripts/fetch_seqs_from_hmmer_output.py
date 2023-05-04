import os
import re
from glob import glob
from collections import defaultdict, Counter
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
@click.option('-h', '--hmmsearch-fn', required=True)
@click.option('-f', '--protein-fasta-dir', required=True)
@click.option('-d', '--output-fasta-dir', required=False, default=None)
@click.option('-1', '--output-one-file', required=False, default=None)
@click.option('-s', '--keep-species-name', default=False, is_flag=True)
def main(hmmsearch_fn, protein_fasta_dir, output_fasta_dir,
         output_one_file, keep_species_name):
    if output_fasta_dir is None and output_one_file is None:
        raise ValueError('must specify either -d or -1')
    fasta_fns = glob_file_type(protein_fasta_dir, '.fa')
    seq_ids_to_fetch = defaultdict(set)
    with open(hmmsearch_fn) as f:
        for record in f:
            record = re.split('\s+', record.strip(), 18)
            seq_id = record[0]
            genus, species, seq_id = seq_id.split('_', 2)
            species_name = f'{genus}_{species}'
            seq_ids_to_fetch[species_name].add(seq_id)

    if output_fasta_dir:
        for species_name, species_seq_ids in seq_ids_to_fetch.items():
            assert species_name in fasta_fns.keys(), species_name
            species_seqs = read_fasta(fasta_fns[species_name])
            output_fn = os.path.join(output_fasta_dir,f'{species_name}.fa')
            with open(output_fn, 'w') as f:
                for seq_id in species_seq_ids:
                    seq = species_seqs[seq_id]
                    f.write(f'>{seq_id}\n{seq}\n')
    else:
        with open(output_one_file, 'w') as f:
            for species_name, species_seq_ids in seq_ids_to_fetch.items():
                assert species_name in fasta_fns.keys(), species_name
                species_seqs = read_fasta(fasta_fns[species_name])
                for seq_id in species_seq_ids:
                    seq = species_seqs[seq_id]
                    if keep_species_name:
                        f.write(f'>{species_name}_{seq_id}\n{seq}\n')
                    else:
                        f.write(f'>{seq_id}\n{seq}\n')


if __name__ == '__main__':
    main()