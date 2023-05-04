

def est_genome_size(wc, input):
    size = 0
    with open(f'{input.genome}.fai') as fai:
        for chrom in fai:
            size += int(chrom.split()[1])
    return int(round(size / 1_000_000))


rule infernal_search:
    input:
        genome=f'../../annotations/{CLADE}/genomic/{{species}}.fa',
        query='../../annotations/cms/{snrna}.cm'
    output:
        tsv='cm_hits/{species}.{snrna}_hits.tsv',
        bed='cm_hits/{species}.{snrna}_hits.bed',
    conda:
        'env_yamls/infernal.yml'
    params:
        genome_size=est_genome_size
    threads: 5
    shell:
        '''
        cmsearch --cpu {threads} --noali -E 0.01 \
          -Z {params.genome_size} \
          --tblout {output.tsv} \
          -o /dev/null \
          {input.query} {input.genome}
        grep  -v "^#" {output.tsv} | \
          awk -v OFS='\t' '{{ \
            if ($10 == "+"){{s=$8; e=$9}} \
            else {{s=$9; e=$8}}; \
            print $1, s, e, "{wildcards.species}|" $16, int($15), $10 \
            }}' | \
        sort -k5,5nr  \
          > {output.bed}
        '''


rule bedtools_getfasta:
    input:
        genome=f'../../annotations/{CLADE}/genomic/{{species}}.fa',
        bed='cm_hits/{species}.{snrna}_hits.bed'
    output:
        fasta='cm_hits/{species}.{snrna}_hits.fasta'
    params:
        n_best=5
    conda:
        'env_yamls/bedtools.yml'
    shell:
        '''
        bedtools getfasta -s -name -bed {input.bed} \
          -fi {input.genome} \
          -fo {output}
        '''


rule sequence_conservation_analysis:
    input:
        u6_aln=expand(
            'cm_hits/{species}.{{snrna}}_hits.fasta',
            species=glob_wildcards(f'../../annotations/{CLADE}/genomic/{{species}}.fa').species
        ),
        tree='orthogroups/species_tree_rooted.txt',
    output:
        svg='figures/{snrna}_aln.svg'
    conda:
        'env_yamls/nb_rpy2_trees.yaml'
    log:
        notebook='notebook_processed/{snrna}_sequence_conservation_analysis.py.ipynb'
    notebook:
        'notebook_templates/{wildcards.snrna}_sequence_conservation_analysis.py.ipynb'

