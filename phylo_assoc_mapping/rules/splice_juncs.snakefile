CLADE = config['clade']


def splice_junc_input(wc):
    if wc.species_name in config['funannotate']['trusted_annotations']:
        return f'../../annotations/{CLADE}/gtf/{wc.species_name}.gtf'
    else:
        return f'annotations/gtf/{wc.species_name}.gtf'
    return prot_fns


def which_script(wc):
    if wc.species_name in config['funannotate']['trusted_annotations']:
        return '../../scripts/gtf_to_cds_intron_bed.py'
    else:
        return '../../scripts/funannotate_gtf_to_cds_intron_bed.py'


rule annot_to_juncs:
    input:
        gtf=splice_junc_input
    output:
        bed='juncs/all_cds_juncs/{species_name}.bed'
    params:
        script=which_script
    conda:
        'env_yamls/pysam.yaml'
    shell:
        '''
        python {params.script} \
          -g {input.gtf} -o {output.bed}
        '''


def splice_junc_filter_input(wc):
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
        'juncs': expand('juncs/all_cds_juncs/{species_name}.bed', species_name=species_names),
        'prot': prot_fns,
    }


rule label_splice_juncs_by_msa_col:
    input:
        unpack(splice_junc_filter_input)
    output:
        juncs=expand(
            'juncs/labelled_cds_juncs/{species_name}.conserved_introns.bed',
            species_name=glob_wildcards(f'../../annotations/{CLADE}/genomic/{{sp}}.fa').sp
        )
    params:
        input_bed_dir='juncs/all_cds_juncs',
        output_bed_dir='juncs/labelled_cds_juncs'
    conda:
        'env_yamls/python_mafft.yaml'
    threads: 36
    shell:
        '''
        mkdir -p $TMPDIR/prot_fastas
        cp {input.prot} $TMPDIR/prot_fastas
        python ../../scripts/conserved_introns.py -p {threads} \
          -o {input.og} \
          -b {params.input_bed_dir} \
          -f $TMPDIR/prot_fastas \
          -d {params.output_bed_dir}
        '''


rule bootstrap_splicing_matrices:
    input:
        fasta=f'../../annotations/{CLADE}/genomic/{{species}}.fa',
        bed='juncs/labelled_cds_juncs/{species}.conserved_introns.bed',
    output:
        pssm='pssms/{species}.{signal_type}.tsv'
    conda:
        'env_yamls/pysam.yaml'
    shell:
        '''
        python ../../scripts/introns_to_motifs.py --nboot 100 \
          -f {input.fasta} -b {input.bed} \
          -o {output.pssm} -s {wildcards.signal_type}
        '''
        
        
rule infer_ancestral_matrices:
    input:
        pssms=expand(
            'pssms/{species}.{{signal_type}}.tsv',
            species=glob_wildcards(f'../../annotations/{CLADE}/genomic/{{species}}.fa').species,
        ),
        tree=ancient('orthogroups/species_tree_rooted.txt'),
    output:
        expand(
            'pssms_inferred/{species}.{{signal_type}}.tsv',
            species=glob_wildcards(f'../../annotations/{CLADE}/genomic/{{species}}.fa').species,
        )
    conda:
        'env_yamls/py_trees.yaml'
    threads: 48
    shell:
        '''
        python ../../scripts/sankoff_motif_analysis.py -n {threads} -i 0.05 \
          -t {input.tree} \
          -s {wildcards.signal_type} \
          -o pssms_inferred \
          {input.pssms}
        '''


rule draw_phenograms:
    input:
        pssms=expand(
            'pssms_inferred/{species}.{{signal_type}}.tsv',
            species=glob_wildcards(f'../../annotations/{CLADE}/genomic/{{species}}.fa').species,
        ),
        tree='orthogroups/species_tree.txt',
        key_ogs='key_prot_ogs/orthogroups/orthogroups_n0.tsv',
        key_prot_ids=ancient(config['key_prot_table_fn']),
    output:
        phenogram='figures/{signal_type}.{phenotype}.phenogram.svg',
        collapsed_tree='figures/{signal_type}.{phenotype}.tree_with_stacked_bars.svg',
        collapsed_tree_cbar='figures/{signal_type}.{phenotype}.tree_with_stacked_bars_cbar.svg',
        full_circular_tree='figures/{signal_type}.{phenotype}.circular_tree.svg'
    conda:
        'env_yamls/nb_rpy2_trees.yaml'
    log:
        notebook='notebook_processed/draw_phenogram.{signal_type}.{phenotype}.py.ipynb'
    notebook:
        'notebook_templates/draw_phenogram.py.ipynb'


rule draw_phenotype_variation:
    input:
        pssms=expand(
            'pssms_inferred/{species}.{{signal_type}}.tsv',
            species=glob_wildcards(f'../../annotations/{CLADE}/genomic/{{species}}.fa').species,
        )
    output:
        pairwise_plot='figures/{signal_type}.pairwise_phenotype_variation.svg',
        ambig_base_plot='figures/{signal_type}.ambig_base_phenotype_variation.svg',
    conda:
        'env_yamls/nb_rpy2_trees.yaml'
    log:
        notebook='notebook_processed/draw_phenotype_variation.{signal_type}.py.ipynb'
    notebook:
        'notebook_templates/draw_phenotype_variation.py.ipynb'


rule og_pgls:
    input:
        pssms=expand(
            'pssms_inferred/{species}.{{signal_type}}.tsv',
            species=glob_wildcards(f'../../annotations/{CLADE}/genomic/{{species}}.fa').species,
        ),
        tree='orthogroups/species_tree_rooted.txt',
        key_ogs='key_prot_ogs/orthogroups/orthogroups_n0.tsv',
        key_prot_ids=ancient(config['key_prot_table_fn']),
    output:
        svg='figures/{signal_type}.{phenotype}.volcano.svg',
        og_phenograms=expand(
            'figures/og_phenograms/{og_id}.{{signal_type}}.{{phenotype}}.phenogram.svg',
            og_id=config['plot_ogs']
        )
    conda:
        'env_yamls/nb_rpy2_trees.yaml'
    log:
        notebook='notebook_processed/og_pgls.{signal_type}.{phenotype}.py.ipynb'
    notebook:
        'notebook_templates/og_pgls.py.ipynb'


rule phmm_pgls:
    input:
        prots=expand(
            'key_prot_ogs/key_og_seqs/{species_name}.fa',
            species_name=glob_wildcards(
                f'../../annotations/{CLADE}/genomic/{{species_name}}.fa'
            ).species_name
        ),
        pssms=expand(
            'pssms_inferred/{species}.{{signal_type}}.tsv',
            species=glob_wildcards(f'../../annotations/{CLADE}/genomic/{{species}}.fa').species,
        ),
        phmm='../../annotations/hmms/{og_id}.hmm',
        tree='orthogroups/species_tree_rooted.txt',
        key_ogs='key_prot_ogs/orthogroups/orthogroups_n0.tsv',
        key_prot_ids=ancient(config['key_prot_table_fn']),
    output:
        stemplot='figures/{og_id}.{signal_type}.{phenotype}.phmm_pgls_stemplot.svg'
    conda:
        'env_yamls/nb_rpy2_trees.yaml'
    log:
        notebook='notebook_processed/phmm_pgls.{og_id}.{signal_type}.{phenotype}.py.ipynb'
    notebook:
        'notebook_templates/phmm_pgls.py.ipynb'


rule u56rho_analysis:
    input:
        juncs=expand(
            'juncs/labelled_cds_juncs/{species}.conserved_introns.bed',
            species=glob_wildcards(f'../../annotations/{CLADE}/genomic/{{species}}.fa').species,
        ),
        fastas=expand(
            f'../../annotations/{CLADE}/genomic/{{species}}.fa',
            species=glob_wildcards(f'../../annotations/{CLADE}/genomic/{{species}}.fa').species,
        ),
        tree='orthogroups/species_tree_rooted.txt',
        key_ogs='key_prot_ogs/orthogroups/orthogroups_n0.tsv',
        key_prot_ids=ancient(config['key_prot_table_fn']),
    output:
        scatter='figures/u56rho.scatter.svg',
        volcano='figures/u56rho.volcano.svg',
    conda:
        'env_yamls/nb_rpy2_trees.yaml'
    log:
        notebook='notebook_processed/u56rho_analysis.py.ipynb'
    notebook:
        'notebook_templates/u56rho_analysis.py.ipynb'