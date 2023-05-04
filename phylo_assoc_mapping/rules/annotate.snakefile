import os
from glob import glob

CLADE = config['clade']


rule funannotate_clean:
    input:
        genome=f'../../annotations/{CLADE}/genomic/{{species_name}}.fa'
    output:
        genome='annotations/genomic/{species_name}.sm.fa'
    conda:
        'env_yamls/funannotate.yaml'
    threads: 12
    shell:
        '''
        funannotate clean -i {input.genome} -o {output.genome}.cleaned.fa
        funannotate sort --simplify -i {output.genome}.cleaned.fa -o {output.genome}.sorted.fa
        funannotate mask --cpus {threads} -i {output.genome}.sorted.fa -o {output}
        rm {output.genome}.cleaned.fa {output.genome}.sorted.fa
        '''


def get_seed_species(wc):
    # don't rerun annotation for trusted assemblies
    assert wc.species_name not in config['funannotate']['trusted_annotations']
    for seed_species, species_names in config['funannotate']['seed_species'].items():
        if wc.species_name in species_names:
            return seed_species
    else:
        raise ValueError(f'Couldn\'t find seed species for {wc.species_name}')


rule funannotate_predict:
    input:
        genome='annotations/genomic/{species_name}.sm.fa'
    output:
        gff3='annotations/gff3/{species_name}.gff3',
        prot='annotations/prot/{species_name}.fa',
    params:
        funannotate_db=os.path.abspath('../../annotations/funannotate'),
        seed_species=get_seed_species,
        ploidy=lambda wc: config['funannotate']['ploidy'].get(wc.species_name, 1),
        intron_len=lambda wc: config['funannotate'].get('intron_len', 3000),
        orgtype=lambda wc: config['funannotate'].get('orgtype', 'fungus')
    conda:
        'env_yamls/funannotate.yaml'
    threads: 12
    shell:
        '''
        GENEMARK_PATH=$(readlink -f ../../scripts/gmes_linux_64_4/) \
          funannotate predict --cpus {threads} \
            --name {wildcards.species_name} \
            --organism {params.orgtype} \
            --ploidy {params.ploidy} \
            --max_intronlen {params.intron_len} \
            -d {params.funannotate_db} \
            --busco_seed_species {params.seed_species} \
            --tmpdir $TMPDIR \
            -i {input.genome} \
            -o $TMPDIR/{wildcards.species_name} \
            -s {wildcards.species_name}
        mv $TMPDIR/{wildcards.species_name}/predict_results/{wildcards.species_name}.gff3 {output.gff3}
        mv $TMPDIR/{wildcards.species_name}/predict_results/{wildcards.species_name}.proteins.fa {output.prot}
        '''


rule gff3_to_gtf:
    input:
        'annotations/gff3/{species_name}.gff3'
    output:
        'annotations/gtf/{species_name}.gtf'
    conda:
        'env_yamls/gffread.yaml'
    threads: 1
    shell:
        '''
        gffread -FCT --no-pseudo {input} > {output}
        '''
