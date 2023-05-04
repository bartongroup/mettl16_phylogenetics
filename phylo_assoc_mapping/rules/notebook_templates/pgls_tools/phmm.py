import shutil
import re
import os
from copy import copy
import subprocess as sp
import tempfile

import pandas as pd


blosum62 = {
    'A': {
        'A': 4.0, 'R': -1.0, 'N': -2.0, 'D': -2.0, 'C': 0.0,
        'Q': -1.0, 'E': -1.0, 'G': 0.0, 'H': -2.0, 'I': -1.0,
        'L': -1.0, 'K': -1.0, 'M': -1.0, 'F': -2.0, 'P': -1.0,
        'S': 1.0, 'T': 0.0, 'W': -3.0, 'Y': -2.0, 'V': 0.0,
        'B': -2.0, 'J': -1.0, 'Z': -1.0, 'X': -1.0, '*': -4.0,
    },
    'R': {
        'A': -1.0, 'R': 5.0, 'N': 0.0, 'D': -2.0, 'C': -3.0,
        'Q': 1.0, 'E': 0.0, 'G': -2.0, 'H': 0.0, 'I': -3.0,
        'L': -2.0, 'K': 2.0, 'M': -1.0, 'F': -3.0, 'P': -2.0,
        'S': -1.0, 'T': -1.0, 'W': -3.0, 'Y': -2.0, 'V': -3.0,
        'B': -1.0, 'J': -2.0, 'Z': 0.0, 'X': -1.0, '*': -4.0,
    },
    'N': {
        'A': -2.0, 'R': 0.0, 'N': 6.0, 'D': 1.0, 'C': -3.0,
        'Q': 0.0, 'E': 0.0, 'G': 0.0, 'H': 1.0, 'I': -3.0,
        'L': -3.0, 'K': 0.0, 'M': -2.0, 'F': -3.0, 'P': -2.0,
        'S': 1.0, 'T': 0.0, 'W': -4.0, 'Y': -2.0, 'V': -3.0,
        'B': 4.0, 'J': -3.0, 'Z': 0.0, 'X': -1.0, '*': -4.0,
    },
    'D': {
        'A': -2.0, 'R': -2.0, 'N': 1.0, 'D': 6.0, 'C': -3.0,
        'Q': 0.0, 'E': 2.0, 'G': -1.0, 'H': -1.0, 'I': -3.0,
        'L': -4.0, 'K': -1.0, 'M': -3.0, 'F': -3.0, 'P': -1.0,
        'S': 0.0, 'T': -1.0, 'W': -4.0, 'Y': -3.0, 'V': -3.0,
        'B': 4.0, 'J': -3.0, 'Z': 1.0, 'X': -1.0, '*': -4.0,
    },
    'C': {
        'A': 0.0, 'R': -3.0, 'N': -3.0, 'D': -3.0, 'C': 9.0,
        'Q': -3.0, 'E': -4.0, 'G': -3.0, 'H': -3.0, 'I': -1.0,
        'L': -1.0, 'K': -3.0, 'M': -1.0, 'F': -2.0, 'P': -3.0,
        'S': -1.0, 'T': -1.0, 'W': -2.0, 'Y': -2.0, 'V': -1.0,
        'B': -3.0, 'J': -1.0, 'Z': -3.0, 'X': -1.0, '*': -4.0,
    },
    'Q': {
        'A': -1.0, 'R': 1.0, 'N': 0.0, 'D': 0.0, 'C': -3.0,
        'Q': 5.0, 'E': 2.0, 'G': -2.0, 'H': 0.0, 'I': -3.0,
        'L': -2.0, 'K': 1.0, 'M': 0.0, 'F': -3.0, 'P': -1.0,
        'S': 0.0, 'T': -1.0, 'W': -2.0, 'Y': -1.0, 'V': -2.0,
        'B': 0.0, 'J': -2.0, 'Z': 4.0, 'X': -1.0, '*': -4.0,
    },
    'E': {
        'A': -1.0, 'R': 0.0, 'N': 0.0, 'D': 2.0, 'C': -4.0,
        'Q': 2.0, 'E': 5.0, 'G': -2.0, 'H': 0.0, 'I': -3.0,
        'L': -3.0, 'K': 1.0, 'M': -2.0, 'F': -3.0, 'P': -1.0,
        'S': 0.0, 'T': -1.0, 'W': -3.0, 'Y': -2.0, 'V': -2.0,
        'B': 1.0, 'J': -3.0, 'Z': 4.0, 'X': -1.0, '*': -4.0,
    },
    'G': {
        'A': 0.0, 'R': -2.0, 'N': 0.0, 'D': -1.0, 'C': -3.0,
        'Q': -2.0, 'E': -2.0, 'G': 6.0, 'H': -2.0, 'I': -4.0,
        'L': -4.0, 'K': -2.0, 'M': -3.0, 'F': -3.0, 'P': -2.0,
        'S': 0.0, 'T': -2.0, 'W': -2.0, 'Y': -3.0, 'V': -3.0,
        'B': -1.0, 'J': -4.0, 'Z': -2.0, 'X': -1.0, '*': -4.0,
    },
    'H': {
        'A': -2.0, 'R': 0.0, 'N': 1.0, 'D': -1.0, 'C': -3.0,
        'Q': 0.0, 'E': 0.0, 'G': -2.0, 'H': 8.0, 'I': -3.0,
        'L': -3.0, 'K': -1.0, 'M': -2.0, 'F': -1.0, 'P': -2.0,
        'S': -1.0, 'T': -2.0, 'W': -2.0, 'Y': 2.0, 'V': -3.0,
        'B': 0.0, 'J': -3.0, 'Z': 0.0, 'X': -1.0, '*': -4.0,
    },
    'I': {
        'A': -1.0, 'R': -3.0, 'N': -3.0, 'D': -3.0, 'C': -1.0,
        'Q': -3.0, 'E': -3.0, 'G': -4.0, 'H': -3.0, 'I': 4.0,
        'L': 2.0, 'K': -3.0, 'M': 1.0, 'F': 0.0, 'P': -3.0,
        'S': -2.0, 'T': -1.0, 'W': -3.0, 'Y': -1.0, 'V': 3.0,
        'B': -3.0, 'J': 3.0, 'Z': -3.0, 'X': -1.0, '*': -4.0,
    },
    'L': {
        'A': -1.0, 'R': -2.0, 'N': -3.0, 'D': -4.0, 'C': -1.0,
        'Q': -2.0, 'E': -3.0, 'G': -4.0, 'H': -3.0, 'I': 2.0,
        'L': 4.0, 'K': -2.0, 'M': 2.0, 'F': 0.0, 'P': -3.0,
        'S': -2.0, 'T': -1.0, 'W': -2.0, 'Y': -1.0, 'V': 1.0,
        'B': -4.0, 'J': 3.0, 'Z': -3.0, 'X': -1.0, '*': -4.0,
    },
    'K': {
        'A': -1.0, 'R': 2.0, 'N': 0.0, 'D': -1.0, 'C': -3.0,
        'Q': 1.0, 'E': 1.0, 'G': -2.0, 'H': -1.0, 'I': -3.0,
        'L': -2.0, 'K': 5.0, 'M': -1.0, 'F': -3.0, 'P': -1.0,
        'S': 0.0, 'T': -1.0, 'W': -3.0, 'Y': -2.0, 'V': -2.0,
        'B': 0.0, 'J': -3.0, 'Z': 1.0, 'X': -1.0, '*': -4.0,
    },
    'M': {
        'A': -1.0, 'R': -1.0, 'N': -2.0, 'D': -3.0, 'C': -1.0,
        'Q': 0.0, 'E': -2.0, 'G': -3.0, 'H': -2.0, 'I': 1.0,
        'L': 2.0, 'K': -1.0, 'M': 5.0, 'F': 0.0, 'P': -2.0,
        'S': -1.0, 'T': -1.0, 'W': -1.0, 'Y': -1.0, 'V': 1.0,
        'B': -3.0, 'J': 2.0, 'Z': -1.0, 'X': -1.0, '*': -4.0,
    },
    'F': {
        'A': -2.0, 'R': -3.0, 'N': -3.0, 'D': -3.0, 'C': -2.0,
        'Q': -3.0, 'E': -3.0, 'G': -3.0, 'H': -1.0, 'I': 0.0,
        'L': 0.0, 'K': -3.0, 'M': 0.0, 'F': 6.0, 'P': -4.0,
        'S': -2.0, 'T': -2.0, 'W': 1.0, 'Y': 3.0, 'V': -1.0,
        'B': -3.0, 'J': 0.0, 'Z': -3.0, 'X': -1.0, '*': -4.0,
    },
    'P': {
        'A': -1.0, 'R': -2.0, 'N': -2.0, 'D': -1.0, 'C': -3.0,
        'Q': -1.0, 'E': -1.0, 'G': -2.0, 'H': -2.0, 'I': -3.0,
        'L': -3.0, 'K': -1.0, 'M': -2.0, 'F': -4.0, 'P': 7.0,
        'S': -1.0, 'T': -1.0, 'W': -4.0, 'Y': -3.0, 'V': -2.0,
        'B': -2.0, 'J': -3.0, 'Z': -1.0, 'X': -1.0, '*': -4.0,
    },
    'S': {
        'A': 1.0, 'R': -1.0, 'N': 1.0, 'D': 0.0, 'C': -1.0,
        'Q': 0.0, 'E': 0.0, 'G': 0.0, 'H': -1.0, 'I': -2.0,
        'L': -2.0, 'K': 0.0, 'M': -1.0, 'F': -2.0, 'P': -1.0,
        'S': 4.0, 'T': 1.0, 'W': -3.0, 'Y': -2.0, 'V': -2.0,
        'B': 0.0, 'J': -2.0, 'Z': 0.0, 'X': -1.0, '*': -4.0,
    },
    'T': {
        'A': 0.0, 'R': -1.0, 'N': 0.0, 'D': -1.0, 'C': -1.0,
        'Q': -1.0, 'E': -1.0, 'G': -2.0, 'H': -2.0, 'I': -1.0,
        'L': -1.0, 'K': -1.0, 'M': -1.0, 'F': -2.0, 'P': -1.0,
        'S': 1.0, 'T': 5.0, 'W': -2.0, 'Y': -2.0, 'V': 0.0,
        'B': -1.0, 'J': -1.0, 'Z': -1.0, 'X': -1.0, '*': -4.0,
    },
    'W': {
        'A': -3.0, 'R': -3.0, 'N': -4.0, 'D': -4.0, 'C': -2.0,
        'Q': -2.0, 'E': -3.0, 'G': -2.0, 'H': -2.0, 'I': -3.0,
        'L': -2.0, 'K': -3.0, 'M': -1.0, 'F': 1.0, 'P': -4.0,
        'S': -3.0, 'T': -2.0, 'W': 11.0, 'Y': 2.0, 'V': -3.0,
        'B': -4.0, 'J': -2.0, 'Z': -2.0, 'X': -1.0, '*': -4.0,
    },
    'Y': {
        'A': -2.0, 'R': -2.0, 'N': -2.0, 'D': -3.0, 'C': -2.0,
        'Q': -1.0, 'E': -2.0, 'G': -3.0, 'H': 2.0, 'I': -1.0,
        'L': -1.0, 'K': -2.0, 'M': -1.0, 'F': 3.0, 'P': -3.0,
        'S': -2.0, 'T': -2.0, 'W': 2.0, 'Y': 7.0, 'V': -1.0,
        'B': -3.0, 'J': -1.0, 'Z': -2.0, 'X': -1.0, '*': -4.0,
    },
    'V': {
        'A': 0.0, 'R': -3.0, 'N': -3.0, 'D': -3.0, 'C': -1.0,
        'Q': -2.0, 'E': -2.0, 'G': -3.0, 'H': -3.0, 'I': 3.0,
        'L': 1.0, 'K': -2.0, 'M': 1.0, 'F': -1.0, 'P': -2.0,
        'S': -2.0, 'T': 0.0, 'W': -3.0, 'Y': -1.0, 'V': 4.0,
        'B': -3.0, 'J': 2.0, 'Z': -2.0, 'X': -1.0, '*': -4.0,
    },
    'B': {
        'A': -2.0, 'R': -1.0, 'N': 4.0, 'D': 4.0, 'C': -3.0,
        'Q': 0.0, 'E': 1.0, 'G': -1.0, 'H': 0.0, 'I': -3.0,
        'L': -4.0, 'K': 0.0, 'M': -3.0, 'F': -3.0, 'P': -2.0,
        'S': 0.0, 'T': -1.0, 'W': -4.0, 'Y': -3.0, 'V': -3.0,
        'B': 4.0, 'J': -3.0, 'Z': 0.0, 'X': -1.0, '*': -4.0,
    },
    'J': {
        'A': -1.0, 'R': -2.0, 'N': -3.0, 'D': -3.0, 'C': -1.0,
        'Q': -2.0, 'E': -3.0, 'G': -4.0, 'H': -3.0, 'I': 3.0,
        'L': 3.0, 'K': -3.0, 'M': 2.0, 'F': 0.0, 'P': -3.0,
        'S': -2.0, 'T': -1.0, 'W': -2.0, 'Y': -1.0, 'V': 2.0,
        'B': -3.0, 'J': 3.0, 'Z': -3.0, 'X': -1.0, '*': -4.0,
    },
    'Z': {
        'A': -1.0, 'R': 0.0, 'N': 0.0, 'D': 1.0, 'C': -3.0,
        'Q': 4.0, 'E': 4.0, 'G': -2.0, 'H': 0.0, 'I': -3.0,
        'L': -3.0, 'K': 1.0, 'M': -1.0, 'F': -3.0, 'P': -1.0,
        'S': 0.0, 'T': -1.0, 'W': -2.0, 'Y': -2.0, 'V': -2.0,
        'B': 0.0, 'J': -3.0, 'Z': 4.0, 'X': -1.0, '*': -4.0,
    },
    'X': {
        'A': -1.0, 'R': -1.0, 'N': -1.0, 'D': -1.0, 'C': -1.0,
        'Q': -1.0, 'E': -1.0, 'G': -1.0, 'H': -1.0, 'I': -1.0,
        'L': -1.0, 'K': -1.0, 'M': -1.0, 'F': -1.0, 'P': -1.0,
        'S': -1.0, 'T': -1.0, 'W': -1.0, 'Y': -1.0, 'V': -1.0,
        'B': -1.0, 'J': -1.0, 'Z': -1.0, 'X': -1.0, '*': -4.0,
    },
    '*': {
        'A': -4.0, 'R': -4.0, 'N': -4.0, 'D': -4.0, 'C': -4.0,
        'Q': -4.0, 'E': -4.0, 'G': -4.0, 'H': -4.0, 'I': -4.0,
        'L': -4.0, 'K': -4.0, 'M': -4.0, 'F': -4.0, 'P': -4.0,
        'S': -4.0, 'T': -4.0, 'W': -4.0, 'Y': -4.0, 'V': -4.0,
        'B': -4.0, 'J': -4.0, 'Z': -4.0, 'X': -4.0, '*': 1.0,
    },
}


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


def sanitise_protein_seqs(seqs):
    cleaned = {
        seq_id: re.sub('[^ACDEFGHIKLMNPQRSTVWY]', 'X', seq)
        for seq_id, seq in seqs.items()
    }
    return cleaned


def check_cmd(cmd, cmd_name):
    if shutil.which(cmd) is None:
        raise RuntimeError(f'{cmd_name} does not appear to be installed')


def hmm_build(hmm_fn, seqs, max_iter=1000, fast_method_threshold=500):
    check_cmd('mafft', 'MAFFT')
    check_cmd('hmmbuild', 'HMMER')
    check_cmd('hmmpress', 'HMMER')
    tmp_fa_fh, tmp_fa_fn = tempfile.mkstemp()
    if len(seqs) <= fast_method_threshold:
        align_cmd = f"mafft --localpair --maxiterate {max_iter} {tmp_fa_fn}"
    else:
        align_cmd = f"mafft --retree 2 --maxiterate {max_iter} {tmp_fa_fn}"
    build_cmd = align_cmd + f' | hmmbuild -n hmm --informat AFA {hmm_fn} -'
    with open(tmp_fa_fn, 'w') as tmp:
        for seq_id, seq in sanitise_protein_seqs(seqs).items():
            tmp.write(f">{seq_id}\n{seq}\n")
    build_proc = sp.run(build_cmd, shell=True, check=True,
                        stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    press_proc = sp.run(
        f'hmmpress -f {hmm_fn}',
        shell=True,
        check=True,
        stdout=sp.DEVNULL,
        stderr=sp.DEVNULL)
    os.close(tmp_fa_fh)
    os.remove(tmp_fa_fn)


def hmm_consensus(hmm_fn):
    check_cmd('hmmemit', 'HMMER')
    proc = sp.run(
        f'hmmemit -c {hmm_fn}',
        shell=True,
        check=True,
        stderr=sp.DEVNULL,
        stdout=sp.PIPE
    )
    seq = _parse_fasta(iter(proc.stdout.decode().split('\n')))
    assert len(seq) == 1
    return list(seq.values())[0]


def _hmmalign(seqs, hmm_fn):
    check_cmd('hmmalign', 'HMMER')
    tmp_fh, tmp_fn = tempfile.mkstemp()
    with open(tmp_fn, 'w') as tmp:
        for seq_id, seq in sanitise_protein_seqs(seqs).items():
            tmp.write(f">{seq_id}\n{seq}\n")
    proc = sp.run(
        f'hmmalign --outformat afa --trim {hmm_fn} {tmp_fn}',
        shell=True,
        check=True,
        stdout=sp.PIPE,
        stderr=sp.DEVNULL
    )
    os.close(tmp_fh)
    os.remove(tmp_fn)
    return _parse_fasta(iter(proc.stdout.decode().split('\n')))


def hmm_align(fasta_fn, hmm_fn=None):
    seqs = read_fasta(fasta_fn)
    if hmm_fn is None:
        hmm_fh, hmm_fn = tempfile.mkstemp()
        hmm_build(hmm_fn, seqs)
        remove_hmm = True
    else:
        remove_hmm = False
    hmm_cons = hmm_consensus(hmm_fn)
    seqs = copy(seqs)
    seqs['hmm_consensus'] = hmm_cons
    aln = _hmmalign(seqs, hmm_fn)
    hmm_cons_aln = aln.pop('hmm_consensus')
    if remove_hmm:
        os.close(hmm_fh)
        os.remove(hmm_fn)
    return build_consensus_table(aln, hmm_cons_aln)


def build_consensus_table(seqs, consensus_seq, method='match'):
    cons_table = []
    aa_table = []
    index = []
    for seq_id, seq in seqs.items():
        genus, species, *_ = seq_id.split('_')
        species = f'{genus}_{species}'
        seq_cons_row = []
        for i, n in enumerate(seq):
            if consensus_seq[i] not in '.-':
                if method == 'match':
                    seq_cons_row.append(int(n == consensus_seq[i]))
                else:
                    if n in '.-':
                        seq_cons_row.append(-4.0)
                    else:
                        seq_cons_row.append(blosum62[n][consensus_seq[i]])
        index.append(species)
        cons_table.append(seq_cons_row)
    columns = [n for n in consensus_seq if n not in '.-']
    columns = [f'{i}_{n}' for i, n in enumerate(columns)]
    cons_table = pd.DataFrame(cons_table, index=index, columns=columns)
    cons_table = cons_table.groupby(level=0).first()
    return cons_table.T


def iter_positions(cons_table, min_species_with_var=5):
    for pos_id, pos in cons_table.iterrows():
        if min_species_with_var <= pos.sum() <= len(pos) - min_species_with_var:
            yield pos_id, pos