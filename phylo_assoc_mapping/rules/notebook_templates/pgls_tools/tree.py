import itertools as it

import numpy as np
import ete3


def get_covariance_from_tree(tree, species):
    root = tree.get_tree_root()
    cov = np.zeros(shape=(len(species), len(species)))

    for i, j in it.product(range(len(species)), repeat=2):
        s1 = species[i]
        s2 = species[j]
        if s1 == s2:
            c = root.get_distance(s1)
        else:
            lca = tree.get_common_ancestor(s1, s2)
            c = root.get_distance(lca)
        cov[i, j] = c
    return cov


def label_tree_pav(node, og_pav):
    if not node.is_leaf():
        assert len(node.children) == 2
        left, right = node.children
        label_tree_pav(left, og_pav)
        label_tree_pav(right, og_pav)
        if left.pav or right.pav:
            node.add_feature('pav', 1)
        else:
            node.add_feature('pav', 0)
    else:
        node.add_feature('pav', int(bool(og_pav[node.name])))


def count_losses(node, og_pav):
    label_tree_pav(node, og_pav)
    losses = []
    for n in node.traverse():
        if n.pav == 0:
            if not n.is_root() and n.up.pav:
                losses.append(n.name)
    return losses


def dollo_gain_loss(tree, og_pav):
    tree = tree.copy()
    # first get the most recent common ancestor of species with an og copy
    species_with_og_copy = og_pav[og_pav == 1].index.tolist()
    if not len(species_with_og_copy):
        # no proteins!
        return None, []
    if len(species_with_og_copy) == 1:
        # only one species has a copy of this protein, special case
        return species_with_og_copy[0], []
    else:
        lca = tree.get_common_ancestor(*species_with_og_copy)
        return lca.name, count_losses(lca, og_pav)


def read_tree(tree_fn):
    tree = ete3.Tree(tree_fn, format=1)
    # name nodes of tree:
    i = 0
    for node in tree.traverse('preorder'):
        if not node.name:
            node.name = f'n{i}'
            i += 1
    return tree

