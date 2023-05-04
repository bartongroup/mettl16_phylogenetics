import os
import itertools as it
from collections import defaultdict

import numpy as np
import pandas as pd
import ete3
import click
from joblib import Parallel, delayed

from matplotlib import pyplot as plt
from matplotlib_logo import logo


def get_cost_matrix(states, cost_func):
    cmat = np.empty(shape=(len(states), len(states)))
    np.fill_diagonal(cmat, 0)
    for i, j in it.combinations(range(len(states)), r=2):
        c_ij = cost_func(states[i], states[j])
        cmat[i, j] = c_ij
        cmat[j, i] = c_ij
    return cmat
                        

def child_cost(child_costs, cost_matrix):
    costs = (cost_matrix + child_costs)
    return costs.min(1), costs.argmin(1)


def _sankoff_up(node, leaf_states, states, cost_matrix):
    if node.is_leaf():
        cost = np.array([
             0.0 if np.allclose(s, leaf_states[node.name]) else np.inf
            for s in states
        ])
        node.add_feature('cost', cost)
    else:
        left, right = node.children
        _sankoff_up(left, leaf_states, states, cost_matrix)
        _sankoff_up(right, leaf_states, states, cost_matrix)
        left_min_cost, left_argmin_cost = child_cost(left.cost, cost_matrix)
        right_min_cost, right_argmin_cost = child_cost(right.cost, cost_matrix)
        node.add_feature('cost', left_min_cost + right_min_cost)
        node.add_feature('left_traceback', left_argmin_cost)
        node.add_feature('right_traceback', right_argmin_cost)


def _sankoff_down(i, node, states):
    node.add_feature('state', states[i])
    node.del_feature('cost')
    if not node.is_leaf():
        left_i = node.left_traceback[i]
        right_i = node.right_traceback[i]
        node.del_feature('left_traceback')
        node.del_feature('right_traceback')
        left, right = node.children
        _sankoff_down(left_i, left, states)
        _sankoff_down(right_i, right, states)


def sankoff(tree, leaf_states, states, cost_func):
    cost_matrix = get_cost_matrix(states, cost_func)
    _tree = tree.copy()
    _sankoff_up(_tree, leaf_states, states, cost_matrix)
    i = np.argmin(_tree.cost)
    min_cost = _tree.cost[i]
    _sankoff_down(i, _tree, states)
    return _tree, min_cost


def round_retain_sum(x, increment=0.1):
    x = x / increment
    N = np.round(x.sum()).astype(int)
    y = x.astype(int)
    M = np.sum(y)
    K = N - M 
    z = x-y 
    if K!=0:
        idx = np.argsort(z)[-K:]
        y[idx] += 1     
    return y * increment


def sm_cost(sm1, sm2):
    return ((sm1 - sm2) ** 2).sum()


def get_possible_sm(increment=0.1):
    # convert to integer
    tot = int(1 / increment)
    vals = np.arange(0, tot + 1)
    for i, j, k, n in it.product(vals, repeat=4):
        if i + j + k + n == tot:
            yield np.array([i, j, k, n]) / tot


def parallel_infer(n, pssms, tree, npos, possible_states, sm_cost, increment):
    boot_pssms = {sp: p.loc[pd.IndexSlice[n, :]].values for sp, p in pssms.items()}
    anc_states = defaultdict(list)
    click.echo(f'starting boostrap {n}')
    for i in range(npos):
        i_leaf_states = {
            species: round_retain_sum(p[i], increment)
            for species, p in boot_pssms.items()
        }
        labelled_tree, _ = sankoff(
            tree, i_leaf_states,
            possible_states,
            cost_func=sm_cost
        )
        for node in labelled_tree.traverse():
            anc_states[node.name].append(node.state)

    anc_states_df = {}
    for node, states in anc_states.items():
        states = pd.DataFrame(states, columns=['A', 'C', 'G', 'U'])
        states.index.name = 'position'
        anc_states_df[node] = states
    click.echo(f'finished boostrap {n}')

    return anc_states_df


def draw_logo(counts,
              alphabet='rna',
              y_format='bits',
              fig_height=2,
              fig_width_per_base=1,
              cmap=None,
              ax=None,
              pad=0.05,
              base_imgs_path=None):
    alphabet = logo.ALPHABETS[alphabet.lower()]
    #counts = calculate_normalised_counts(seqs, alphabet)
    entropy = logo.calculate_entropy(counts)
    if y_format == 'bits':
        heights = logo.calculate_bits(counts, entropy)
        ylim = heights.sum(1).max()
    elif y_format == 'probability':
        heights = counts
        ylim = 1
    else:
        raise ValueError(
            'y_format "{}" is not recognised, use "bits" or "probability"'.format(y_format))
    if ax is None:
        fig, ax = plt.subplots(figsize=(fig_width_per_base * len(entropy), fig_height))
    imgs = logo._get_base_imgs(logo._get_unambiguous(alphabet), base_imgs_path, cmap)
    for i, position_heights in enumerate(heights):
        order = np.argsort(position_heights)
        base_order = [list(alphabet)[i] for i in order]
        position_heights_sorted = position_heights[order]
        y = 0
        for base_idx, h in zip(base_order, position_heights_sorted):
            if h:
                ax.imshow(imgs[base_idx],
                          origin='upper',
                          extent=[i + pad / 2, i + 1 - pad / 2, y, y + h],
                          interpolation='bilinear')
            y += h
    ax.set_xlim(0, i + 1)
    ax.set_ylim(0, ylim)
    ax.set_ylabel(y_format.capitalize())
    ax.set_xticks([])
    ax.set_aspect('auto')
    return fig, ax


def create_itol_stacked_bar_dataset(states, npos, signal_type, output_dir):
    
    for pos in range(npos):
        itol_dataset = (
f'''DATASET_MULTIBAR
SEPARATOR SPACE
DATASET_LABEL {signal_type}_pos{pos+1}
COLOR #0072b2
FIELD_COLORS #009e73 #0072b2 #f0e442 #d55e00
FIELD_LABELS A C G U
SHOW_INTERNAL 0

LEGEND_TITLE 5SS+4pos
LEGEND_POSITION_X 100
LEGEND_POSITION_Y 100
LEGEND_HORIZONTAL 0
LEGEND_SHAPES 1 1 1 1
LEGEND_COLORS #009e73 #0072b2 #f0e442 #d55e00
LEGEND_LABELS A C G U

DATA
'''
        )
        for node_name, pssm in states.items():
            data = pssm.values[pos]
            itol_dataset += f'{node_name} ' + ' '.join([f'{d:.2f}' for d in data]) + '\n'
        with open(os.path.join(output_dir, f'{signal_type}_pos{pos + 1}.itol_multibar.txt'), 'w') as f:
            f.write(itol_dataset)


@click.command()
@click.argument('pssm-fns', required=True, nargs=-1)
@click.option('-t', '--species-tree-fn', required=True)
@click.option('-o', '--output-dir', required=True)
@click.option('-s', '--signal-type', type=click.Choice(['5ss', '3ss', 'bp']), default='5ss')
@click.option('-i', '--increment', default=0.05)
@click.option('-n', '--n-jobs', default=1)
def intron_seq_pssm(pssm_fns, species_tree_fn, output_dir, signal_type, increment, n_jobs):
    tree = ete3.Tree(species_tree_fn, format=1)
    pssms = {
        os.path.split(fn)[1].split('.')[0]: pd.read_csv(fn, sep='\t', index_col=[0, 1])
        for fn in pssm_fns
    }
    nboots = list(pssms.values())[0].index.get_level_values('bootstrap').max() + 1
    npos = list(pssms.values())[0].index.get_level_values('position').max() + 1
    possible_states = list(get_possible_sm(increment))
    with Parallel(n_jobs) as pool:
        boot_anc_states_list = pool(
            delayed(parallel_infer)(n, pssms, tree, npos, possible_states, sm_cost, increment)
            for n in range(nboots)
        )
    boot_anc_states = defaultdict(dict)
    for n, anc_states in enumerate(boot_anc_states_list):
        for node, states in anc_states.items():
            boot_anc_states[node][n] = states
    boot_anc_states = {sp: pd.concat(states, axis=0, names=['bootstrap'])
                       for sp, states in boot_anc_states.items()}
    av_anc_states = {sp: states.groupby('position').mean() for sp, states in boot_anc_states.items()}
    create_itol_stacked_bar_dataset(av_anc_states, npos, signal_type, output_dir)
    for node, states in boot_anc_states.items():
        output_fn = os.path.join(output_dir, f'{node}.{signal_type}.tsv')
        states.to_csv(output_fn, sep='\t')
        fig, ax = draw_logo(av_anc_states[node].values, y_format='probability')
        plt.savefig(output_fn[:-4] + '.svg')
        plt.close()

if __name__ == '__main__':
    intron_seq_pssm()