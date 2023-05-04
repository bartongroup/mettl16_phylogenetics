CLADE = config['clade']


def orthofinder_input(wc):
    species_names = glob_wildcards(f'../../annotations/{CLADE}/genomic/{{species_name}}.fa').species_name
    prot_fns = []
    for sp in species_names:
        if sp in config['funannotate']['trusted_annotations']:
            prot_fns.append(f'../../annotations/{CLADE}/prot_repr/{sp}.fa')
        else:
            prot_fns.append(f'annotations/prot/{sp}.fa')
    return prot_fns


rule orthofinder:
    input:
        proteomes=orthofinder_input
    output:
        n0='orthogroups/orthogroups_n0.tsv',
        tree='orthogroups/species_tree_rooted.txt',
        tree_with_confidences='orthogroups/species_tree_confidences.txt',
        stats='orthogroups/orthofinder_stats_overall.tsv',
        stats_per_species='orthogroups/orthofinder_stats_per_species.tsv'
    threads: 48
    conda:
        'env_yamls/orthofinder.yaml'
    shell:
        '''
        mkdir -p $TMPDIR/orthofinder_input
        cp {input.proteomes} $TMPDIR/orthofinder_input
        USE_MEM=1 ../../scripts/OrthoFinder/orthofinder.py -t 44 -a 4 \
          -f $TMPDIR/orthofinder_input \
          -n orthofinder \
          -o $TMPDIR/orthofinder_output
        mkdir -p orthogroups
        mv $TMPDIR/orthofinder_output/Results_orthofinder/Phylogenetic_Hierarchical_Orthogroups/N0.tsv {output.n0}
        mv $TMPDIR/orthofinder_output/Results_orthofinder/Species_Tree/SpeciesTree_rooted_node_labels.txt {output.tree}
        mv $TMPDIR/orthofinder_output/Results_orthofinder/Species_Tree/SpeciesTree_rooted.txt {output.tree_with_confidences}
        mv $TMPDIR/orthofinder_output/Results_orthofinder/Comparative_Genomics_Statistics/Statistics_Overall.tsv {output.stats}
        mv $TMPDIR/orthofinder_output/Results_orthofinder/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv {output.stats_per_species}
        '''


rule merge_trees:
    input:
        tree_with_node_labels='orthogroups/species_tree_rooted.txt',
        tree_with_confidences='orthogroups/species_tree_confidences.txt',
    output:
        'orthogroups/species_tree.txt',
    conda:
        'env_yamls/py_trees.yaml'
    shell:
        '''
        python ../../scripts/reformat_tree.py \
          {input.tree_with_node_labels} \
          {input.tree_with_confidences} \
          {output}
        '''