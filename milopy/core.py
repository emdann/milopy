from audioop import add
import logging
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
import re

# rpy2 setup to run edgeR functions
from rpy2.robjects.packages import PackageNotInstalledError, importr
import rpy2.robjects.numpy2ri
from rpy2.robjects import pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import STAP


def make_nhoods(
        adata: AnnData,
        neighbors_key: str = None,
        prop=0.1,
        seed: int = 42):
    '''
    This function randomly samples vertcies on a graph to define neighbourhoods.
    These are then refined by computing the median profile for the neighbourhood
    in reduced dimensional space and selecting the nearest vertex to this
    position. Thus, multiple neighbourhoods may be collapsed down together to
    prevent over-sampling the graph space.

    Params:
    -------
    - adata: AnnData object. Should contain a knn graph in `adata.obsp`
    - neighbors_key: string indicating the key in `adata.obsp` to use as KNN graph. If not specified, 
    `make_nhoods` looks .obsp[‘connectivities’] for connectivities (default storage places for
    `scanpy.pp.neighbors`). If specified, it looks .obsp[.uns[neighbors_key][‘connectivities_key’]] for
    connectivities.
    - prop: fraction of cells to sample for neighbourhood index search (default: 0.1)
    - seed: random seed for cell sampling (default: 42)

    Returns:
    --------
    None, adds in place:
    - `adata.obsm['nhoods']`: a binary matrix of cell to neighbourhood assignments
    - `adata.obs['nhood_ixs_refined']`: a boolean indicating whether a cell is an index for a neighbourhood
    - `adata.obs['kth_distance']`: the distance to the kth nearest neighbour for each index cell (used for SpatialFDR correction)
    - `adata.uns["nhood_neighbors_key"]`: stores the KNN graph key in `adata.obsp` used for neighbourhood construction
    '''
    # Get reduced dim used for KNN graph
    if neighbors_key is None:
        try:
            use_rep = adata.uns["neighbors"]["params"]["use_rep"]
        except KeyError:
            logging.warning('Using X_pca as default embedding')
            use_rep = "X_pca"
        try:
            knn_graph = adata.obsp["connectivities"].copy()
        except KeyError:
            raise KeyError(
                'No "connectivities" slot in adata.obsp -- please run scanpy.pp.neighbors(adata) first'
            )
    else:
        try:
            use_rep = adata.uns["neighbors"]["params"]["use_rep"]
        except KeyError:
            logging.warning('Using X_pca as default embedding')
            use_rep = "X_pca"
        knn_graph = adata.obsp[neighbors_key + "_connectivities"].copy()

    # Get reduced dim
    X_dimred = adata.obsm[use_rep]

    # Sample size
    n_ixs = int(np.round(adata.n_obs * prop))

    # Binarize
    knn_graph[knn_graph != 0] = 1

    #  Sample random vertices
    random.seed(seed)
    random_vertices = random.sample(range(adata.n_obs), k=n_ixs)
    random_vertices.sort()

    ixs_nn = knn_graph[random_vertices, :]

    # Refine sampling
    non_zero_rows = ixs_nn.nonzero()[0]
    non_zero_cols = ixs_nn.nonzero()[1]

    nh_reduced_dims = np.empty(shape=[0, X_dimred.shape[1]])
    refined_vertices = np.empty(shape=[len(random_vertices), ])

    from sklearn.metrics.pairwise import euclidean_distances
    for i in range(len(random_vertices)):
        nh_pos = np.median(
            X_dimred[non_zero_cols[non_zero_rows == i], :], 0).reshape(-1, 1)
        nn_ixs = non_zero_cols[non_zero_rows == i]
        # Find closest real point (amongst nearest neighbors)
        dists = euclidean_distances(
            X_dimred[non_zero_cols[non_zero_rows == i], :], nh_pos.T)
        # Update vertex index
        refined_vertices[i] = nn_ixs[dists.argmin()]

    refined_vertices = np.unique(refined_vertices.astype("int"))
    refined_vertices.sort()

    nhoods = knn_graph[:, refined_vertices]
    adata.obsm['nhoods'] = nhoods

    # Add ixs to adata
    adata.obs["nhood_ixs_random"] = adata.obs_names.isin(
        adata.obs_names[random_vertices])
    adata.obs["nhood_ixs_refined"] = adata.obs_names.isin(
        adata.obs_names[refined_vertices])
    adata.obs["nhood_ixs_refined"] = adata.obs["nhood_ixs_refined"].astype(
        "int")
    adata.obs["nhood_ixs_random"] = adata.obs["nhood_ixs_random"].astype("int")
    # Store info on neighbor_key used
    adata.uns["nhood_neighbors_key"] = neighbors_key
    # Store distance to K-th nearest neighbor (used for spatial FDR correction)
    if neighbors_key is None:
        k = adata.uns["neighbors"]["params"]["n_neighbors"]
        knn_dists = adata.obsp["distances"]
    else:
        k = adata.uns[neighbors_key]["params"]["n_neighbors"]
        knn_dists = adata.obsp[neighbors_key + "_distances"]
    nhood_ixs = adata.obs["nhood_ixs_refined"] == 1
    dist_mat = knn_dists[nhood_ixs, :]
    k_distances = dist_mat.max(1).toarray().ravel()
    adata.obs["nhood_kth_distance"] = 0
    adata.obs.loc[adata.obs["nhood_ixs_refined"]
                  == 1, "nhood_kth_distance"] = k_distances


def count_nhoods(
    adata: AnnData,
    sample_col: str,
):
    '''
    Builds a sample-level AnnData object storing the matrix of cell counts per sample per neighbourhood.

    Params:
    -------
    - adata: AnnData object with neighbourhoods defined in `adata.obsm['nhoods']`
    - sample_col: string, column in adata.obs that contains sample information 
    (what should be in the columns of the cell count matrix)

    Returns: 
    --------
    MuData object storing the original (i.e. cell-level) AnnData in `mudata['cells']`
    and the sample-level anndata storing the neighbourhood cell counts in `mudata['samples']`.
    Here:
    - `mudata['samples'].obs_names` are samples (defined from `adata.obs['sample_col']`)
    - `mudata['samples'].var_names` are neighbourhoods 
    - `mudata['samples'].X` is the matrix counting the number of cells from each
    sample in each neighbourhood
    '''
    try:
        nhoods = adata.obsm["nhoods"]
    except KeyError:
        raise KeyError(
            'Cannot find "nhoods" slot in adata.obsm -- please run milopy.make_nhoods(adata)'
        )
    #  Make nhood abundance matrix
    sample_dummies = pd.get_dummies(adata.obs[sample_col])
    all_samples = sample_dummies.columns
    sample_dummies = scipy.sparse.csr_matrix(sample_dummies.values)
    nhood_count_mat = adata.obsm["nhoods"].T.dot(sample_dummies)
    sample_obs = pd.DataFrame(index=all_samples)
    sample_adata = anndata.AnnData(X=nhood_count_mat.T, obs=sample_obs)
    sample_adata.uns["sample_col"] = sample_col
    # Save nhood index info
    sample_adata.var["index_cell"] = adata.obs_names[adata.obs["nhood_ixs_refined"] == 1]
    sample_adata.var["kth_distance"] = adata.obs.loc[adata.obs["nhood_ixs_refined"]
                                                     == 1, "nhood_kth_distance"].values
    # adata.uns["sample_adata"] = sample_adata
    # Make MuData object
    milo_mdata = MuData({"cells": adata, "samples": sample_adata})
    return(milo_mdata)


def DA_nhoods(milo_mdata: MuData,
              design: str,
              model_contrasts: str = None,
              subset_samples: List[str] = None,
              add_intercept: bool = True):
    '''
    Performs differential abundance testing on neighbourhoods using QLF test implementation from edgeR
    (using R code under the hood)

    Params:
    -------
    - milo_mdata: MuData object, output of `count_nhoods`
    - design: formula for the test, following glm syntax from R (e.g. '~ condition'). Terms should be columns in `milo_mdata['cells'].obs`.
    - model_contrasts: A string vector that defines the contrasts used to perform DA testing, following glm syntax from R (e.g. "conditionDisease - conditionControl").
        If no contrast is specified (default), then the last categorical level in condition of interest is used as the test group.
    - subset_samples: subset of samples (obs in `milo_mdata['samples']`) to use for the test
    - add_intercept: whether to include an intercept in the model. If False, this is equivalent to adding + 0 in the design formula. 
    When model_contrasts is specified, this is set to False by default. 

    Returns:
    --------
    None, modifies `milo_mdata['samples']` in place, adding the results of the DA test to `.var`:
    - `logFC` stores the log fold change in cell abundance (coefficient from the GLM)
    - `PValue` stores the p-value for the QLF test before multiple testing correction
    - `SpatialFDR` stores the the p-value adjusted for multiple testing to limit the false discovery rate, 
        calculated with weighted Benjamini-Hochberg procedure
    '''

    # Get data
    try:
        sample_adata = milo_mdata['samples']
    except KeyError:
        raise ValueError(
            "milo_mdata should be a MuData object with two slots: 'cells' and 'samples' - please run milopy.count_nhoods(adata) first")
    adata = milo_mdata['cells']

    # Set up rpy2 to run edgeR
    rpy2.robjects.numpy2ri.activate()
    rpy2.robjects.pandas2ri.activate()
    edgeR = _try_import_bioc_library("edgeR")
    limma = _try_import_bioc_library("limma")
    stats = importr("stats")
    base = importr("base")

    covariates = [x.strip(" ") for x in set(
        re.split('\\+|\\*', design.lstrip("~ ")))]

    # Add covariates used for testing to sample_adata.var
    sample_col = sample_adata.uns["sample_col"]
    try:
        sample_obs = adata.obs[covariates + [sample_col]].drop_duplicates()
    except KeyError:
        missing_cov = [
            x for x in covariates if x not in sample_adata.obs.columns]
        raise KeyError(
            'Covariates {c} are not columns in adata.obs'.format(
                c=" ".join(missing_cov))
        )

    sample_obs = sample_obs[covariates + [sample_col]]
    sample_obs.index = sample_obs[sample_col].astype("str")

    try:
        assert sample_obs.loc[sample_adata.obs_names].shape[0] == len(
            sample_adata.obs_names)
    except:
        raise ValueError(
            "Covariates cannot be unambiguously assigned to each sample -- each sample value should match a single covariate value")
    sample_adata.obs = sample_obs.loc[sample_adata.obs_names]
    # Get design dataframe
    try:
        design_df = sample_adata.obs[covariates]
    except KeyError:
        missing_cov = [
            x for x in covariates if x not in sample_adata.obs.columns]
        raise KeyError(
            'Covariates {c} are not columns in adata.uns["sample_adata"].obs'.format(
                c=" ".join(missing_cov))
        )

    # Get count matrix
    count_mat = sample_adata.X.T.toarray()
    lib_size = count_mat.sum(0)

    # Filter out samples with zero counts
    keep_smp = lib_size > 0

    # Subset samples
    if subset_samples is not None:
        keep_smp = keep_smp & sample_adata.obs_names.isin(subset_samples)
        design_df = design_df[keep_smp]
        for i, e in enumerate(design_df.columns):
            if design_df.dtypes[i].name == 'category':
                design_df[e] = design_df[e].cat.remove_unused_categories()

    # Filter out nhoods with zero counts
    # (they can appear after sample filtering)
    keep_nhoods = count_mat[:, keep_smp].sum(1) > 0

    # Define model matrix
    if not add_intercept or model_contrasts is not None:
        design = design + ' + 0'
    model = stats.model_matrix(object=stats.formula(
        design), data=design_df)

    # Fit NB-GLM
    dge = edgeR.DGEList(
        counts=count_mat[keep_nhoods, :][:, keep_smp], lib_size=lib_size[keep_smp])
    dge = edgeR.calcNormFactors(dge, method="TMM")
    dge = edgeR.estimateDisp(dge, model)
    fit = edgeR.glmQLFit(dge, model, robust=True)

    # Test
    n_coef = model.shape[1]
    if model_contrasts is not None:
        r_str = '''
        get_model_cols <- function(design_df, design){
            m = model.matrix(object=formula(design), data=design_df)
            return(colnames(m))
        }
        '''
        get_model_cols = STAP(r_str, "get_model_cols")
        model_mat_cols = get_model_cols.get_model_cols(design_df, design)
        model_df = pd.DataFrame(model)
        model_df.columns = model_mat_cols
        try:
            mod_contrast = limma.makeContrasts(
                contrasts=model_contrasts, levels=model_df)
        except:
            raise ValueError(
                "Model contrasts must be in the form 'A-B' or 'A+B'")
        res = base.as_data_frame(edgeR.topTags(edgeR.glmQLFTest(
            fit, contrast=mod_contrast), sort_by='none', n=np.inf))
    else:
        res = base.as_data_frame(edgeR.topTags(
            edgeR.glmQLFTest(fit, coef=n_coef), sort_by='none', n=np.inf))
    res = rpy2.robjects.conversion.rpy2py(res)
    if not isinstance(res, pd.DataFrame):
        res = pd.DataFrame(res)

    # Save outputs
    res.index = sample_adata.var_names[keep_nhoods]
    if any([x in sample_adata.var.columns for x in res.columns]):
        sample_adata.var = sample_adata.var.drop(res.columns, axis=1)
    sample_adata.var = pd.concat([sample_adata.var, res], axis=1)

    # Run Graph spatial FDR correction
    _graph_spatialFDR(
        sample_adata, neighbors_key=adata.uns["nhood_neighbors_key"])


def _graph_spatialFDR(sample_adata, neighbors_key=None):
    '''
    FDR correction weighted on inverse of connectivity of neighbourhoods.
    The distance to the k-th nearest neighbor is used as a measure of connectivity.
    '''

    # use 1/connectivity as the weighting for the weighted BH adjustment from Cydar
    w = 1/sample_adata.var['kth_distance']
    w[np.isinf(w)] = 0

    # Computing a density-weighted q-value.
    pvalues = sample_adata.var["PValue"]
    keep_nhoods = ~pvalues.isna()  #  Filtering in case of test on subset of nhoods
    o = pvalues[keep_nhoods].argsort()
    pvalues = pvalues[keep_nhoods][o]
    w = w[keep_nhoods][o]

    adjp = np.zeros(shape=len(o))
    adjp[o] = (sum(w)*pvalues/np.cumsum(w))[::-1].cummin()[::-1]
    adjp = np.array([x if x < 1 else 1 for x in adjp])

    ## Store in anndata
    sample_adata.var["SpatialFDR"] = np.nan
    sample_adata.var.loc[keep_nhoods, "SpatialFDR"] = adjp

## -- UTILS -- ##


def _try_import_bioc_library(name):
    try:
        _r_lib = importr(name)
        return _r_lib
    except PackageNotInstalledError:
        raise RuntimeError(
            f"Install Bioconductor library `{name!r}` first as `BiocManager::install({name!r}).`"
        )
