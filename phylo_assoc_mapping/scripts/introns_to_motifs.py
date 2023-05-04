import os
import sys
import re
from collections import defaultdict, Counter
import itertools as it

import numpy as np
import pandas as pd

import pysam
import click


RC = str.maketrans('ACGU', 'UGCA')


IUPAC = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'U': {'U'},
    'S': {'G', 'C'},
    'W': {'A', 'U'},
    'R': {'A', 'G'},
    'Y': {'C', 'U'},
    'M': {'A', 'C'},
    'K': {'G', 'U'},
    'B': {'C', 'G', 'U'},
    'D': {'A', 'G', 'U'},
    'H': {'A', 'C', 'U'},
    'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'U'}
}

IUPAC_REV = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'U': 'U',
    'CG': 'S',
    'AU': 'W',
    'AG': 'R',
    'CU': 'Y',
    'AC': 'M',
    'GU': 'K',
    'CGU': 'B',
    'AGU': 'D',
    'ACU': 'H',
    'ACG': 'V',
    'ACGU': 'N'
}


def rev_comp(seq):
    return seq.translate(RC)[::-1]


def predict_branchpoint(seq):
    bp_motifs = list(re.finditer('[CU]U[AG]A[CU]', seq))
    if not len(bp_motifs):
        return None
    else:
        bp_motif = bp_motifs[-1]
        start = bp_motif.start()
        lpad = abs(min(0, start - 3))
        start = max(0, start - 3)
        end = bp_motif.end() + 2
        rpad = max(end - len(seq), 0)
        bp = seq[start: end]
        bp = lpad * 'N' + bp + rpad * 'N'
        return bp


def read_intron_seqs(bed_fn, fasta_fn, signal_type):
    introns = set()
    with open(bed_fn) as bed:
        for record in bed:
            chrom, start, end, _, _, strand, *_ = record.strip().split('\t')
            start, end = int(start), int(end)
            introns.add((chrom, start, end, strand))
    seqs = []
    with pysam.FastaFile(fasta_fn) as fasta:
        for chrom, start, end, strand in introns:
            if signal_type == '5ss':
                left_pos = start - 3 if strand == '+' else end - 7
                right_pos = start + 7 if strand == '+' else end + 3
                seq = fasta.fetch(chrom, left_pos, right_pos).upper().replace('T', 'U')
                if strand == '-':
                    seq = rev_comp(seq)
            elif signal_type == '3ss':
                left_pos = end - 5 if strand == '+' else start - 3
                right_pos = end + 3 if strand == '+' else start + 5
                seq = fasta.fetch(chrom, left_pos, right_pos).upper().replace('T', 'U')
                if strand == '-':
                    seq = rev_comp(seq)
            elif signal_type == 'bp':
                intron_seq = fasta.fetch(chrom, start, end).upper().replace('T', 'U')
                if strand == '-':
                    intron_seq = rev_comp(intron_seq)
                seq = predict_branchpoint(intron_seq)
            if seq is not None:
                seqs.append(seq)
    return seqs


def calculate_normalised_counts(seqs):
    seqlen = len(seqs[0])
    counts = np.zeros((seqlen, 4), dtype='f')
    alphabet =  {
        'A': np.array([1, 0, 0, 0]),
        'C': np.array([0, 1, 0, 0]),
        'G': np.array([0, 0, 1, 0]),
        'U': np.array([0, 0, 0, 1]),
        'R': np.array([0.5, 0., 0.5, 0.]),
        'Y': np.array([0., 0.5, 0., 0.5]),
        'S': np.array([0., 0.5, 0.5, 0.]),
        'W': np.array([0.5, 0., 0., 0.5]),
        'K': np.array([0., 0., 0.5, 0.5]),
        'M': np.array([0.5, 0.5, 0., 0.]),
        'B': np.array([0., 1./3, 1./3, 1./3]),
        'D': np.array([1./3, 0., 1./3, 1./3]),
        'H': np.array([1./3, 1./3, 0., 1./3]),
        'V': np.array([1./3, 1./3, 1./3, 0.]),
        'N': np.array([0.25, 0.25, 0.25, 0.25]),
    }
    for seq in seqs:
        if len(seq) != seqlen:
            raise ValueError('Not all sequences are of the same length')
        for i, base in enumerate(seq):
            try:
                counts[i] += alphabet[base]
            except KeyError:
                continue
    pssm = pd.DataFrame(counts / counts.sum(1)[:, np.newaxis], columns=['A', 'C', 'G', 'U'])
    pssm.index.name = 'position'
    return pssm


@click.command()
@click.option('-f', '--fasta-fn', required=True)
@click.option('-b', '--bed-fn', required=True)
@click.option('-o', '--output-fn', required=True)
@click.option('-s', '--signal-type', type=click.Choice(['5ss', '3ss', 'bp']), default='5ss')
@click.option('--nboot', type=int, default=25)
def intron_seq_pssm(fasta_fn, bed_fn, output_fn, signal_type, nboot):
    seqs = read_intron_seqs(bed_fn, fasta_fn, signal_type)
    pssms = {}
    for n in range(nboot):
        idx = np.random.randint(0, len(seqs), len(seqs))
        seqs_sample = [seqs[i] for i in idx]
        pssms[n] = calculate_normalised_counts(seqs_sample)
    pssms = pd.concat(pssms, axis=0, names=['bootstrap'])
    pssms.to_csv(output_fn, sep='\t')


if __name__ == '__main__':
    intron_seq_pssm()
