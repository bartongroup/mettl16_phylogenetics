CLADE = config['clade']


def build_hmm_input(wc):
    species_names = glob_wildcards(
        f'../../annotations/{CLADE}/genomic/{{species_name}}.fa'
    ).species_name
    prot_fns = []
    for sp in species_names:
        if sp in config['funannotate']['trusted_annotations']:
            prot_fns.append(f'../../annotations/{CLADE}/prot_repr/{sp}.fa')
        else:
            prot_fns.append(f'annotations/prot/{sp}.fa')
    return {
        'og': 'orthogroups/orthogroups_n0.tsv',
        'key_prot_ids': ancient(config['key_prot_table_fn']),
        'prot': prot_fns,
    }


checkpoint build_hmms:
    input:
        unpack(build_hmm_input)
    output:
        score_thresh='key_prot_ogs/hmms/hmm_scoring_thresholds.tsv'
    params:
        output_dir='key_prot_ogs/hmms/'
    conda:
        'env_yamls/python_mafft.yaml'
    threads: 24
    shell:
        '''
        mkdir -p $TMPDIR/prot_fastas
        cp {input.prot} $TMPDIR/prot_fastas
        python ../../scripts/create_hmms.py -p {threads} \
          -o {input.og} \
          -k {input.key_prot_ids} \
          -f $TMPDIR/prot_fastas \
          -d {params.output_dir}
        '''


rule sixframe_translate:
    input:
        genome=f'../../annotations/{CLADE}/genomic/{{species_name}}.fa'
    output:
        sixframe='annotations/sixframe/{species_name}.fa'
    params:
        min_length=30
    conda:
        'env_yamls/python_mafft.yaml'
    shell:
        '''
        esl-translate -l {params.min_length} {input.genome} > {output.sixframe}
        '''


def search_6frame_input(wc):
    species_names = glob_wildcards(f'../../annotations/{CLADE}/genomic/{{species_name}}.fa').species_name
    prot_fns = []
    for sp in species_names:
        if sp in config['funannotate']['trusted_annotations']:
            prot_fns.append(f'../../annotations/{CLADE}/prot_repr/{sp}.fa')
        else:
            prot_fns.append(f'annotations/sixframe/{sp}.fa')
    return prot_fns


def get_score_threshold(wc):
    hmm_score_table_fn = checkpoints.build_hmms.get().output.score_thresh
    with open(hmm_score_table_fn) as f:
        #skip header
        next(f)
        for record in f:
            _, hmm_name, score_threshold = record.split()
            if hmm_name == wc.hmm_name:
                return float(score_threshold)


rule search_6frame:
    input:
        sixframe=search_6frame_input
    output:
        hmm_hits='key_prot_ogs/hmm_hits/{hmm_name}.tsv',
    params:
        score_threshold=get_score_threshold
    conda:
        'env_yamls/python_mafft.yaml'
    threads: 24
    shell:
        '''
        for FN in {input.sixframe}; do
          BASENAME="${{FN##*/}}"
          SP="${{BASENAME%.fa}}"
          sed -e "s/>/>${{SP}}_/" $FN
        done |
        hmmsearch --cpu {threads} -o /dev/null --tblout /dev/stdout \
          -T {params.score_threshold} \
          key_prot_ogs/hmms/{wildcards.hmm_name}.hmm - | \
        grep -v "^#" > {output.hmm_hits}
        '''


def agg_6frame_input(wc):
    hmm_score_table_fn = checkpoints.build_hmms.get().output.score_thresh
    hmm_names = []
    with open(hmm_score_table_fn) as f:
        #skip header
        next(f)
        for record in f:
            _, hmm_name, score_threshold = record.split()
            hmm_names.append(hmm_name)

    return expand(
        'key_prot_ogs/hmm_hits/{hmm_name}.tsv',
        hmm_name=hmm_names
    )
        

rule aggregate_6frame_results:
    input:
        hmm_hits=agg_6frame_input
    output:
        'key_prot_ogs/all_hmm_hits.tsv'
    shell:
        '''
        cat {input.hmm_hits} | grep -v "^#" > {output}
        '''


def extract_seqs_input(wc):
    species_names = glob_wildcards(
        f'../../annotations/{CLADE}/genomic/{{species_name}}.fa'
    ).species_name
    prot_fns = []
    for sp in species_names:
        if sp in config['funannotate']['trusted_annotations']:
            prot_fns.append(f'../../annotations/{CLADE}/prot_repr/{sp}.fa')
        else:
            prot_fns.append(f'annotations/sixframe/{sp}.fa')
    return {
        'hmm_hits': 'key_prot_ogs/all_hmm_hits.tsv',
        'prot': prot_fns,
    }


rule extract_seqs:
    input:
        unpack(extract_seqs_input)
    output:
        expand(
            'key_prot_ogs/key_og_seqs/{species_name}.fa',
            species_name=glob_wildcards(
                f'../../annotations/{CLADE}/genomic/{{species_name}}.fa'
            ).species_name
        )
    conda:
        'env_yamls/python_mafft.yaml'
    shell:
        '''
        mkdir -p $TMPDIR/prot_fastas
        cp {input.prot} $TMPDIR/prot_fastas
        python ../../scripts/fetch_seqs_from_hmmer_output.py \
          -h {input.hmm_hits} \
          -f $TMPDIR/prot_fastas \
          -d key_prot_ogs/key_og_seqs
        '''


rule reorthofinder_key_ogs:
    input:
        prots=expand(
            'key_prot_ogs/key_og_seqs/{species_name}.fa',
            species_name=glob_wildcards(
                f'../../annotations/{CLADE}/genomic/{{species_name}}.fa'
            ).species_name
        ),
        tree='orthogroups/species_tree_rooted.txt',
    output:
        n0='key_prot_ogs/orthogroups/orthogroups_n0.tsv',
        gene_trees=directory('key_prot_ogs/orthogroups/gene_trees'),
    threads: 36
    conda:
        'env_yamls/orthofinder.yaml'
    shell:
        '''
        mkdir -p $TMPDIR/orthofinder_input
        cp {input.prots} $TMPDIR/orthofinder_input
        USE_MEM=1 ../../scripts/OrthoFinder/orthofinder.py -t 22 -a 2 \
          -s {input.tree} \
          -f $TMPDIR/orthofinder_input \
          -n orthofinder \
          -o $TMPDIR/orthofinder_output
        mkdir -p orthogroups
        mv $TMPDIR/orthofinder_output/Results_orthofinder/Phylogenetic_Hierarchical_Orthogroups/N0.tsv {output.n0}
        mv $TMPDIR/orthofinder_output/Results_orthofinder/Resolved_Gene_Trees {output.gene_trees}
        '''
