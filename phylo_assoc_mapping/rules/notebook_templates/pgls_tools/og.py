import warnings

import numpy as np
import pandas as pd

from .tree import dollo_gain_loss


def get_key_ogs(key_prot_ids, og):
    og_ids = []
    for _, gene_symbol, prot_id, species_name in key_prot_ids.itertuples():
        try:
            og_id = og.query(f'{species_name}.str.contains("{prot_id}")').index[0]
        except IndexError:
            warnings.warn(f'Could not find orthogroup for "{gene_symbol}"')
            og_id = np.nan
        og_ids.append(og_id)
    key_prot_ids['og_id'] = og_ids
    og_idx = []
    new_og_ids = []
    for og_id, group in key_prot_ids.groupby('og_id'):
        og_idx.append(og_id)
        new_og_ids.append('_'.join(group.gene_symbol))
    og_key = og.loc[og_idx]
    og_key.index = new_og_ids
    return og_key


def read_orthofinder_n0(n0_fn, key_prot_fn):
    og = pd.read_csv(n0_fn, sep='\t', index_col='HOG')
    og = og.iloc[:, 2:].fillna('')
    key_prot_table = pd.read_csv(
        key_prot_fn, sep='\t', names=['gene_symbol', 'prot_id', 'species']
    )
    og = get_key_ogs(key_prot_table, og)
    og_cnv = og.applymap(lambda g: len(g.split(', ')) if len(g) else 0)
    return og, og_cnv


def iter_og(n0_fn, key_prot_fn, tree, species,
            allow_gain=False, min_losses=2, min_prot_per_og=5):
    og, og_cnv = read_orthofinder_n0(n0_fn, key_prot_fn)
    og_cnv = og_cnv[species]
    og_pav = (og_cnv != 0).astype(int)
    root_name = tree.get_tree_root().name
    cluster_id = 0
    for _, group in og_pav.groupby(by=species):
        group_repr = group.iloc[0].copy()
        og_ids = group.index.tolist()
        if group_repr.sum() > min_prot_per_og:
            gain, losses = dollo_gain_loss(tree, group_repr)
            if allow_gain or gain == root_name:
                if len(losses) > min_losses:
                    cluster_name = f'OGCLUSTER{cluster_id}'
                    group_repr = group_repr.rename(cluster_name)
                    yield (cluster_name, og_ids, group_repr,
                           int(gain != root_name), len(losses), group_repr.sum())
                    cluster_id += 1