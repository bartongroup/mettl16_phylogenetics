import numpy as  np
import pandas as pd
from scipy import stats
from joblib import Parallel, delayed
from statsmodels.stats.multitest import fdrcorrection

from .og import iter_og
from .phmm import hmm_align, iter_positions
from .tree import read_tree, get_covariance_from_tree
from .pheno import pssms_to_pheno, load_bootpheno
from .pgls import bootstrap_phylogenetic_gls


def _og_parallel_pgls(pheno, cov, exog_data):
    og_cluster_id, og_ids, og, gain_type, nlosses, nprot = exog_data
    coef, coef_lower, coef_upper, p_val = bootstrap_phylogenetic_gls(og, pheno, cov)
    return [
        og_cluster_id, ';'.join(og_ids),
        gain_type, nlosses, nprot,
        coef[0], coef_lower[0], coef_upper[0], p_val[0]
    ]


def pssm_og_pgls(n0_fn, key_prot_fn, pssm_fns, tree_fn, pos, nom_bases, denom_bases,
                 allow_gain=False, min_losses=2, min_prot_per_og=4, njobs=1):
    # load tree and generate covariance matrix
    tree = read_tree(tree_fn)
    species = [n.name for n in tree.get_leaves()]
    cov = get_covariance_from_tree(tree, species)

    # load phenotype information
    pheno = pssms_to_pheno(pssm_fns, pos, nom_bases, denom_bases)
    pheno = pheno.loc[:, species]

    # create iterator of orthogroup clusters (clusters contain OGs with same presense/absence variation)
    og_iter = iter_og(
        n0_fn, key_prot_fn, tree, species,
        allow_gain=allow_gain,
        min_losses=min_losses,
        min_prot_per_og=min_prot_per_og
    )
    res = Parallel(njobs)(
        delayed(_og_parallel_pgls)(pheno, cov, og_data) for og_data in og_iter
    )    
    res = pd.DataFrame(
        res, columns=['cluster_id', 'og_ids', 'gain',
                      'nlosses', 'nprot',
                      'coef', 'coef_lower', 'coef_upper',
                      'p_val']
    )
    _, res['fdr'], *_ = fdrcorrection(res.p_val)
    return res


def bootpheno_og_pgls(n0_fn, key_prot_fn, bootpheno_fn, tree_fn,
                      allow_gain=False, min_losses=2, min_prot_per_og=4, njobs=1):
    # load tree and generate covariance matrix
    tree = read_tree(tree_fn)
    species = [n.name for n in tree.get_leaves()]
    cov = get_covariance_from_tree(tree, species)

    # load phenotype information
    pheno = load_bootpheno(bootpheno_fn)
    pheno = pheno.loc[:, species]

    # create iterator of orthogroup clusters (clusters contain OGs with same presense/absence variation)
    og_iter = iter_og(
        n0_fn, key_prot_fn, tree, species,
        allow_gain=allow_gain,
        min_losses=min_losses,
        min_prot_per_og=min_prot_per_og
    )
    res = Parallel(njobs)(
        delayed(_og_parallel_pgls)(pheno, cov, og_data) for og_data in og_iter
    )    
    res = pd.DataFrame(
        res, columns=['cluster_id', 'og_ids', 'gain',
                      'nlosses', 'nprot',
                      'coef', 'coef_lower', 'coef_upper',
                      'p_val']
    )
    _, res['fdr'], *_ = fdrcorrection(res.p_val)
    return res


def _phmm_parallel_pgls(pheno, cov, exog_data):
    pos_id, pos_var = exog_data
    coef, coef_lower, coef_upper, p_val = bootstrap_phylogenetic_gls(pos_var, pheno, cov)
    return [
        pos_id,
        coef[0], coef_lower[0], coef_upper[0], p_val[0]
    ]


def pssm_phmm_pgls(fasta_fn, phmm_fn, pssm_fns, tree_fn, pos, nom_bases, denom_bases,
        min_species_with_var=5, njobs=1):
    # load tree and generate covariance matrix
    pos_var = hmm_align(fasta_fn, phmm_fn)
    species = pos_var.columns.values
    tree = read_tree(tree_fn)
    cov = get_covariance_from_tree(tree, species)

    # load phenotype information
    pheno = pssms_to_pheno(pssm_fns, pos, nom_bases, denom_bases)
    pheno = pheno.loc[:, species]

    res = Parallel(njobs)(
        delayed(_phmm_parallel_pgls)(pheno, cov, pos_data)
        for pos_data in iter_positions(pos_var, min_species_with_var)
    )    
    res = pd.DataFrame(
        res, columns=['pos_id', 'coef',
                      'coef_lower', 'coef_upper',
                      'p_val']
    )
    _, res['fdr'], *_ = fdrcorrection(res.p_val)
    return res