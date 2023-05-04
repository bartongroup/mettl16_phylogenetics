import re
from operator import itemgetter
import warnings

import numpy as np
import click


def get_gtf_attribute(gtf_record, attribute):
    try:
        attr = re.search(f'(?:^|\s){attribute} "(.+?)";', gtf_record[8]).group(1)
    except AttributeError:
        raise ValueError(
            f'Could not parse attribute {attribute} '
            f'from GTF with feature type {gtf_record[2]}'
        )
    return attr


def gtf_cds_parser(gtf_fn):
    gtf_records = {}
    with open(gtf_fn) as gtf:
        for i, record in enumerate(gtf):
            if record.startswith('#') or record.startswith('track'):
                continue
            record = record.split('\t')
            chrom, _, feat_type, start, end, _, strand= record[:7]
            start = int(start) - 1
            end = int(end)
            if feat_type == 'CDS':
                gene_id = get_gtf_attribute(record, 'gene_id')
                transcript_id = get_gtf_attribute(record, 'transcript_id')
                idx = (chrom, gene_id, transcript_id, strand)
                if idx not in gtf_records:
                    gtf_records[idx] = []
                gtf_records[idx].append([start, end])
    for (chrom, gene_id, transcript_id, strand), invs in gtf_records.items():
        invs.sort(key=itemgetter(0))
        cds = np.array(invs)
        yield chrom, gene_id, transcript_id, strand, cds

        
def exons_to_introns(exons):
    return np.stack([exons[:-1, 1], exons[1:, 0]], axis=1)


def calculate_intron_prot_pos(exon_sizes):
    exon_csizes = np.cumsum(exon_sizes)[:-1]
    prot_pos, prot_frame = np.divmod(exon_csizes, 3)
    return prot_pos, prot_frame


def get_intron_protein_pos(cds, strand):
    exon_sizes = cds[:, 1] - cds[:, 0]
    if strand == '+':
        prot_pos, prot_frame = calculate_intron_prot_pos(exon_sizes)
    else:
        prot_pos, prot_frame = calculate_intron_prot_pos(exon_sizes[::-1])
        prot_pos = prot_pos[::-1]
        prot_frame = prot_frame[::-1]
    return prot_pos, prot_frame


@click.command()
@click.option('-g', '--gtf-fn', required=True)
@click.option('-o', '--output-bed-fn', required=True)
def gtf_to_protein_introns(gtf_fn, output_bed_fn):
    with open(output_bed_fn, 'w') as bed:
        for chrom, gene_id, transcript_id, strand, cds in gtf_cds_parser(gtf_fn):
            introns = exons_to_introns(cds)
            prot_pos, prot_frame = get_intron_protein_pos(cds, strand)
            for (start, end), ppos, pframe in zip(introns, prot_pos, prot_frame):
                bed_record =(
                    f'{chrom}\t{start:d}\t{end:d}\t'
                    f'{transcript_id}|{ppos:d}|{pframe:d}\t'
                    f'.\t{strand}\n'
                )
                bed.write(bed_record)


if __name__ == '__main__':
    gtf_to_protein_introns()