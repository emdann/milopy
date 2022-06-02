import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse
import anndata
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple
from anndata import AnnData
from scipy.sparse import csr_matrix
import random

## -- NHOOD EXPRESSION -- ##


def add_nhood_expression(
    adata: AnnData,
    layer: str = None,
):
    '''
    Calculates the mean expression in neighbourhoods of each feature in `adata.X` or
    `adata.layers[layer]` (if layer is not None).

    Params:
    -------
    - adata: AnnData object
    - layer: which data layer to use as expression matrix (default: None, uses `adata.X`)

    Returns:
    -------
    Updates adata in place to store the matrix of average expression in each neighbourhood in `adata.uns["nhood_adata"].obsm['expr']`
    '''
    try:
        nhood_adata = adata.uns["nhood_adata"]
    except KeyError:
        raise KeyError(
            'Cannot find "nhood_adata" slot in adata.uns -- please run milopy.make_nhoods_adata(adata)'
        )

    # Get gene expression matrix
    if layer is None:
        X = adata.X
        expr_id = "expr"
    else:
        X = adata.layers[layer]
        expr_id = "expr_" + layer

    # Aggregate over nhoods -- taking the mean
    nhoods_X = X.T.dot(adata.obsm["nhoods"])
    nhoods_X = csr_matrix(nhoods_X / adata.obsm["nhoods"].toarray().sum(0))
    adata.uns["nhood_adata"].obsm[expr_id] = nhoods_X.T


## -- NHOOD GRAPH -- ##

def build_nhood_graph(adata, basis="X_umap"):
    '''
    Build graph of neighbourhoods used for visualization of DA results

    Params:
    -------
    - adata: AnnData object
    - basis: string indicating the name of the obsm basis to use to use for layout of neighbourhoods (key in `adata.obsm`)
    '''
    # Add embedding positions
    adata.uns["nhood_adata"].obsm["X_milo_graph"] = adata[adata.obs["nhood_ixs_refined"] == 1].obsm[basis]
    # Add nhood size
    adata.uns["nhood_adata"].obs["Nhood_size"] = np.array(
        adata.obsm["nhoods"].sum(0)).flatten()
    # Add adjacency graph
    adata.uns["nhood_adata"].obsp["nhood_connectivities"] = adata.obsm["nhoods"].T.dot(
        adata.obsm["nhoods"])
    adata.uns["nhood_adata"].uns["nhood"] = {
        "connectivities_key": "nhood_connectivities", "distances_key": ""}

## -- UTILS --- ##


def add_covariate_to_nhoods_var(adata, new_covariates):
    '''
    Add covariate from adata.obs to adata.uns["nhood_adata"].var
    '''
    # Add covariates used for testing to nhood_adata.var
    try:
        nhood_adata = adata.uns["nhood_adata"].copy()
    except KeyError:
        raise KeyError(
            'Cannot find "nhood_adata" slot in adata.uns -- please run milopy.make_nhoods_adata(adata)'
        )

    sample_col = nhood_adata.uns["sample_col"]
    covariates = list(set(
        nhood_adata.var.columns[nhood_adata.var.columns != sample_col].tolist() + new_covariates))
    try:
        nhoods_var = adata.obs[covariates + [sample_col]].drop_duplicates()
    except KeyError:
        missing_cov = [
            x for x in covariates if x not in nhood_adata.var.columns]
        raise KeyError(
            'Covariates {c} are not columns in adata.obs'.format(
                c=" ".join(missing_cov))
        )
    nhoods_var = nhoods_var[covariates + [sample_col]].astype("str")
    nhoods_var.index = nhoods_var[sample_col]
    try:
        assert nhoods_var.loc[nhood_adata.var_names].shape[0] == len(
            nhood_adata.var_names)
    except:
        raise ValueError(
            "Covariates cannot be unambiguously assigned to each sample -- each sample value should match a single covariate value")
    nhood_adata.var = nhoods_var.loc[nhood_adata.var_names]
    adata.uns["nhood_adata"] = nhood_adata


def annotate_nhoods(adata, anno_col):
    '''
    Assigns a categorical label to neighbourhoods, based on the most frequent label
    among cells in each neighbourhood. This can be useful to stratify DA testing
    results by cell types or samples.

    Params:
    -------
    - adata: AnnData object with adata.uns["nhood_adata"]
    - anno_col: string indicating column in adata.obs containing the cell annotations to use for nhood labelling

    Returns:
    --------
    None. Adds in place:
    - `adata.uns["nhood_adata"].obs["nhood_annotation"]`: assigning a label to each nhood
    - `adata.uns["nhood_adata"].obs["nhood_annotation_frac"]` stores the fraciton of cells in the neighbourhood with the assigned label
    - `adata.uns["nhood_adata"].obsm['frac_annotation']`: stores the fraction of cells from each label in each nhood
    - `adata.uns["nhood_adata"].uns["annotation_labels"]`: stores the column names for `adata.uns["nhood_adata"].obsm['frac_annotation']`
    '''
    anno_dummies = pd.get_dummies(adata.obs[anno_col])
    anno_count = adata.obsm["nhoods"].T.dot(
        scipy.sparse.csr_matrix(anno_dummies.values))
    anno_frac = np.array(anno_count/anno_count.sum(1))

    anno_frac = pd.DataFrame(anno_frac,
                             columns=anno_dummies.columns,
                             index=adata.uns["nhood_adata"].obs_names
                             )
    adata.uns["nhood_adata"].obsm["frac_annotation"] = anno_frac.values
    adata.uns["nhood_adata"].uns["annotation_labels"] = anno_frac.columns
    adata.uns["nhood_adata"].obs["nhood_annotation"] = anno_frac.idxmax(1)
    adata.uns["nhood_adata"].obs["nhood_annotation_frac"] = anno_frac.max(1)


def annotate_nhoods_continuous(adata, anno_col):
    '''
    Assigns a continuous value to neighbourhoods, based on mean cell level covariate stored in adata.obs. 
    This can be useful to correlate DA log-foldChanges with continuous covariates such as pseudotime, gene expression scores etc...

    Params:
    -------
    - adata: AnnData object with adata.uns["nhood_adata"]
    - anno_col: string indicating column in adata.obs containing the cell annotations to use for nhood labelling

    Returns:
    --------
    None. Adds in place:
    - `adata.uns["nhood_adata"].obs["nhood_{anno_col}"]`: assigning a continuous value to each nhood
    '''
    anno_val = adata.obsm["nhoods"].T.dot(
        scipy.sparse.csr_matrix(adata.obs[anno_col]).T)

    mean_anno_val = anno_val.toarray()/np.array(adata.obsm["nhoods"].T.sum(1))

    adata.uns["nhood_adata"].obsm[f"nhood_{anno_col}"] = mean_anno_val.values


## -- CHECKS -- ##

def _check_milo_anndata(adata):
    '''
    Checks that anndata object contains the required slots
    '''
    has_nhoods = 'nhoods' in adata.obsm.keys()
    has_nhood_adata = 'nhood_adata' in adata.uns.keys()
