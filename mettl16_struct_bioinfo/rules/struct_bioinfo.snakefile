import csv


checkpoint domain_segmentation:
    input:
        af_pdb=lambda wc: expand('af_structures/{sp}.pdb', sp=config['comps'][wc.comp].values()),
        af_pae=lambda wc: expand('af_structures/{sp}.json', sp=config['comps'][wc.comp].values()),
    output:
        doms='domains/{comp}.doms.csv',
    conda:
        'env_yamls/nb_struct_bioinfo.yaml'
    log:
        notebook='notebook_processed/domain_segmentation.{comp}.py.ipynb'
    notebook:
        'notebook_templates/domain_segmentation.py.ipynb'


def get_dom_fns(wc):
    dom_csv_fn = checkpoints.domain_segmentation.get(comp=wc.comp).output.doms
    input_ = {
        'doms_csv': dom_csv_fn,
        'doms': [],
        'ref_doms': expand('xray_structures/{fn}.pdb',
                           fn=config['reference_xray_structures'].values())
    }
    with open(dom_csv_fn) as f:
        for row in csv.DictReader(f):
            input_['doms'].append(row['dom_fn_highconf'])
    return input_


rule draw_tmscore_heatmap:
    input:
        unpack(get_dom_fns)
    output:
        doms_labelled='domains/{comp}.doms_labelled.csv',
        heatmap='figures/{comp}.tmscore_heatmap.svg',
        boxplot='figures/{comp}.tmscore_boxplot.svg',
    conda:
        'env_yamls/nb_struct_bioinfo.yaml'
    log:
        notebook='notebook_processed/draw_tmscore_heatmap.{comp}.py.ipynb'
    notebook:
        'notebook_templates/draw_tmscore_heatmap.py.ipynb'


rule draw_domain_struct_tree:
    input:
        doms_labelled='domains/{comp}.doms_labelled.csv',
        tree='trees/{comp}.txt'
    output:
        tree='figures/{comp}.domain_struct_tree.svg',
    conda:
        'env_yamls/nb_struct_bioinfo.yaml'
    log:
        notebook='notebook_processed/draw_domain_struct_tree.{comp}.py.ipynb'
    notebook:
        'notebook_templates/draw_domain_struct_tree.py.ipynb'