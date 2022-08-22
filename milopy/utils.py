import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse
import anndata
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple
from anndata import AnnData
from mudata import MuData
from scipy.sparse import csr_matrix
import random

## -- NHOOD EXPRESSION -- ##


def add_nhood_expression(
    milo_mdata: MuData,
    layer: str = None,
):
    '''
    Calculates the mean expression in neighbourhoods of each feature in `milo_mdata['cells'].X` or
    `milo_mdata['cells'][layer]` (if layer is not None).

    Params:
    -------
    - milo_mdata: MuData object
    - layer: which data layer to use as expression matrix (default: None, uses `milo_mdata['cells'].X`)

    Returns:
    -------
    Updates adata in place to store the matrix of average expression in each neighbourhood in `milo_mdata['samples'].varm['expr']`
    '''
    try:
        sample_adata = milo_mdata['samples']
    except KeyError:
        raise ValueError(
            "milo_mdata should be a MuData object with two slots: 'cells' and 'samples' - please run milopy.count_nhoods(adata) first")
    adata = milo_mdata['cells']

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
    sample_adata.varm[expr_id] = nhoods_X.T


## -- NHOOD GRAPH -- ##

def build_nhood_graph(milo_mdata: MuData,
                      basis: str = "X_umap"):
    '''
    Build graph of neighbourhoods used for visualization of DA results

    Params:
    -------
    - milo_mdata: MuData object
    - basis: string indicating the name of the obsm basis to use to use for layout of neighbourhoods (key in `adata.obsm`)
    '''
    adata = milo_mdata['cells']
    # # Add embedding positions
    milo_mdata['samples'].varm["X_milo_graph"] = adata[adata.obs["nhood_ixs_refined"] == 1].obsm[basis]
    # Add nhood size
    milo_mdata['samples'].var["Nhood_size"] = np.array(
        adata.obsm["nhoods"].sum(0)).flatten()
    # Add adjacency graph
    milo_mdata['samples'].varp["nhood_connectivities"] = adata.obsm["nhoods"].T.dot(
        adata.obsm["nhoods"])
    milo_mdata['samples'].varp['nhood_connectivities'].setdiag(0)
    milo_mdata['samples'].varp['nhood_connectivities'].eliminate_zeros()
    milo_mdata['samples'].uns["nhood"] = {
        "connectivities_key": "nhood_connectivities", "distances_key": ""}

## -- UTILS --- ##


def add_covariate_to_nhoods_var(
        milo_mdata: MuData,
        new_covariates: List[str]):
    '''
    Add covariate from adata.obs to adata.uns["nhood_adata"].var
    '''
    try:
        sample_adata = milo_mdata['samples']
    except KeyError:
        raise ValueError(
            "milo_mdata should be a MuData object with two slots: 'cells' and 'samples' - please run milopy.count_nhoods(adata) first")
    adata = milo_mdata['cells']

    sample_col = sample_adata.uns["sample_col"]
    covariates = list(set(
        sample_adata.obs.columns[sample_adata.obs.columns != sample_col].tolist() + new_covariates))
    try:
        sample_obs = adata.obs[covariates + [sample_col]].drop_duplicates()
    except KeyError:
        missing_cov = [
            x for x in covariates if x not in sample_adata.obs.columns]
        raise KeyError(
            'Covariates {c} are not columns in adata.obs'.format(
                c=" ".join(missing_cov))
        )
    sample_obs = sample_obs[covariates + [sample_col]].astype("str")
    sample_obs.index = sample_obs[sample_col]
    try:
        assert sample_obs.loc[sample_adata.obs_names].shape[0] == len(
            sample_adata.obs_names)
    except:
        raise ValueError(
            "Covariates cannot be unambiguously assigned to each sample -- each sample value should match a single covariate value")
    sample_adata.obs = sample_obs.loc[sample_adata.obs_names]
    # adata.uns["nhood_adata"] = nhood_adata


def annotate_nhoods(milo_mdata: MuData,
                    anno_col: str):
    '''
    Assigns a categorical label to neighbourhoods, based on the most frequent label
    among cells in each neighbourhood. This can be useful to stratify DA testing
    results by cell types or samples.

    Params:
    -------
    - milo_mdata: MuData object 
    - anno_col: string indicating column in adata.obs containing the cell annotations to use for nhood labelling

    Returns:
    --------
    None. Adds in place:
    - `milo_mdata['samples'].var["nhood_annotation"]`: assigning a label to each nhood
    - `milo_mdata['samples'].var["nhood_annotation_frac"]` stores the fraciton of cells in the neighbourhood with the assigned label
    - `milo_mdata['samples'].varm['frac_annotation']`: stores the fraction of cells from each label in each nhood
    - `milo_mdata['samples'].uns["annotation_labels"]`: stores the column names for `milo_mdata['samples'].varm['frac_annotation']`
    '''
    try:
        sample_adata = milo_mdata['samples']
    except KeyError:
        raise ValueError(
            "milo_mdata should be a MuData object with two slots: 'cells' and 'samples' - please run milopy.count_nhoods(adata) first")
    adata = milo_mdata['cells']

    # Check value is not numeric
    if pd.api.types.is_numeric_dtype(adata.obs[anno_col]):
        raise ValueError(
            'adata.obs[anno_col] is not of categorical type - please use milopy.utils.annotate_nhoods_continuous for continuous variables')

    anno_dummies = pd.get_dummies(adata.obs[anno_col])
    anno_count = adata.obsm["nhoods"].T.dot(
        scipy.sparse.csr_matrix(anno_dummies.values))
    anno_frac = np.array(anno_count/anno_count.sum(1))

    anno_frac = pd.DataFrame(anno_frac,
                             columns=anno_dummies.columns,
                             index=milo_mdata['samples'].var_names
                             )
    milo_mdata['samples'].varm["frac_annotation"] = anno_frac.values
    milo_mdata['samples'].uns["annotation_labels"] = anno_frac.columns
    milo_mdata['samples'].uns["annotation_obs"] = anno_col
    milo_mdata['samples'].var["nhood_annotation"] = anno_frac.idxmax(1)
    milo_mdata['samples'].var["nhood_annotation_frac"] = anno_frac.max(1)


def annotate_nhoods_continuous(
        milo_mdata: MuData,
        anno_col: str):
    '''
    Assigns a continuous value to neighbourhoods, based on mean cell level covariate stored in adata.obs. 
    This can be useful to correlate DA log-foldChanges with continuous covariates such as pseudotime, gene expression scores etc...

    Params:
    -------
    - adata: MuData object
    - anno_col: string indicating column in adata.obs containing the cell annotations to use for nhood labelling

    Returns:
    --------
    None. Adds in place:
    - `milo_mdata['samples'].var["nhood_{anno_col}"]`: assigning a continuous value to each nhood
    '''
    try:
        sample_adata = milo_mdata['samples']
    except KeyError:
        raise ValueError(
            "milo_mdata should be a MuData object with two slots: 'cells' and 'samples' - please run milopy.count_nhoods(adata) first")
    adata = milo_mdata['cells']

    # Check value is not categorical
    if not pd.api.types.is_numeric_dtype(adata.obs[anno_col]):
        raise ValueError(
            'adata.obs[anno_col] is not of continuous type - please use milopy.utils.annotate_nhoods for categorical variables')

    anno_val = adata.obsm["nhoods"].T.dot(
        scipy.sparse.csr_matrix(adata.obs[anno_col]).T)

    mean_anno_val = anno_val.toarray()/np.array(adata.obsm["nhoods"].T.sum(1))

    milo_mdata['samples'].var[f"nhood_{anno_col}"] = mean_anno_val


# ## -- I/O -- ##
# def write_milo_adata(adata: AnnData,
#                      filepath: str,
#                      **kwargs):
#     '''
#     Save anndata objects after Milo analysis

#     Params:
#     -----
#     - adata: AnnData object with adata.uns["nhood_adata"]
#     - filepath: path to h5ad file to save
#     - **kwargs: arguments passed to scanpy.write_h5ad

#     Returns:
#     -------
#     None, saves 2 AnnData objects in h5ad format. The cell x gene AnnData is saved in filepath.
#     The nhood x sample AnnData is saved in a separate object (location is stored in adata.uns['nhood_adata_filepath'])
#     '''
#     nhood_filepath = filepath.split('.h5ad')[0] + ".nhood_adata.h5ad"
#     adata.uns['nhood_adata_filepath'] = nhood_filepath
#     try:
#         if 'annotation_labels' in adata.uns['nhood_adata'].uns.keys():
#             adata.uns['nhood_adata'].uns['annotation_labels'] = adata.uns['nhood_adata'].uns['annotation_labels'].tolist()
#     except KeyError:
#         raise KeyError(
#             'Cannot find "nhood_adata" slot in adata.uns -- please run milopy.make_nhoods_adata(adata)')
#     nhood_adata = adata.uns["nhood_adata"].copy()
#     nhood_adata.write_h5ad(nhood_filepath, **kwargs)
#     del adata.uns["nhood_adata"]
#     adata.write_h5ad(filepath, **kwargs)


# def read_milo_adata(
#         filepath: str,
#         **kwargs) -> AnnData:
#     '''
#     Read AnnData objects stored after Milo analysis

#     Params:
#     ------
#     - filepath: path to h5ad file storing cell x gene AnnData object
#     - **kwargs: additional arguments passed to scanpy.read_h5ad

#     Returns:
#     -------
#     - AnnData object storing milo slots (adata.obsm['nhoods'], adata.uns['nhood_adata'])
#     '''
#     adata = sc.read_h5ad(filepath, **kwargs)
#     try:
#         nhood_filepath = adata.uns['nhood_adata_filepath']
#     except:
#         raise KeyError('No nhood_adata_file associated to adata')
#     adata.uns["nhood_adata"] = sc.read_h5ad(nhood_filepath, **kwargs)
#     return(adata)
