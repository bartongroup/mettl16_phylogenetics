import os
import numpy as np
import pandas as pd


def pssm_ratio_phenotype(pssm, pos, nom_bases, denom_bases, log=False, pseudo=1e-3):
    nom = pssm.loc[pd.IndexSlice[:, pos], nom_bases].sum(1)
    denom = pssm.loc[pd.IndexSlice[:, pos], denom_bases].sum(1)
    if log:
        nom = np.log2(nom + pseudo)
        denom = np.log2(denom + pseudo)
    return nom - denom


def pssms_to_pheno(pssm_fns, pos, nom_bases, denom_bases, log=False):
    pssms = {
        os.path.split(fn)[1].split('.')[0]: pd.read_csv(fn, sep='\t', index_col=[0, 1])
        for fn in pssm_fns
    }
    pheno = pd.DataFrame.from_dict(
        {sp: pssm_ratio_phenotype(pssm, pos, nom_bases, denom_bases).values
         for sp, pssm in pssms.items()}
    )
    pheno.index.name = 'bootstrap'
    return pheno


def load_bootpheno(bootpheno_fn):
    pheno = pd.read_csv(
        bootpheno_fn, sep='\t', index_col=0
    )
    pheno.index.name = 'bootstrap'
    return pheno