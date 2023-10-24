from audioop import add
import logging
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse
import anndata
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple
from anndata import AnnData
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

    - adata: AnnData object. Should contain a knn graph in `adata.obsp`
    - neighbors_key: string indicating the key in `adata.obsp` to use as KNN graph. If not specified, 
    `make_nhoods` looks .obsp[‘connectivities’] for connectivities (default storage places for
    `scanpy.pp.neighbors`). If specified, it looks .obsp[.uns[neighbors_key][‘connectivities_key’]] for
    connectivities.
    - prop: fraction of cells to sample for neighbourhood index search (default: 0.1)
    - seed: random seed for cell sampling (default: 42)
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
            use_rep = adata.uns[neighbors_key]["params"]["use_rep"]
        except KeyError:
            logging.warning('Using X_pca as default embedding')
            use_rep = "X_pca"
        knn_graph = adata.obsp[neighbors_key + "_connectivities"].copy()

    # Get reduced dim
    if use_rep == 'X':
        X_dimred = adata.X
        if scipy.sparse.issparse(X_dimred)
            X_dimred = X_dimred.A
    else:
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
    - adata
    - sample_col: string, column in adata.obs that contains sample information 
    (what should be in the columns of the nhoodCount matrix)

    Returns: None
    Updated adata.uns slot to contain adata.uns["nhood_adata"], where:
    - adata.uns["nhood_adata"].obs_names are neighbourhoods
    - adata.uns["nhood_adata"].var_names are samples
    - adata.uns["nhood_adata"].X is the matrix counting the number of cells from each
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
    nhood_var = pd.DataFrame(index=all_samples)
    nhood_adata = anndata.AnnData(X=nhood_count_mat, var=nhood_var)
    nhood_adata.uns["sample_col"] = sample_col
    # Save nhood index info
    nhood_adata.obs["index_cell"] = adata.obs_names[adata.obs["nhood_ixs_refined"] == 1]
    nhood_adata.obs["kth_distance"] = adata.obs.loc[adata.obs["nhood_ixs_refined"]
                                                    == 1, "nhood_kth_distance"].values
    adata.uns["nhood_adata"] = nhood_adata


def DA_nhoods(adata, design, model_contrasts=None, subset_samples=None, add_intercept=True):
    '''
    This will perform differential neighbourhood abundance testing (using edgeR under the hood)
    - adata
    - design: formula (terms should be columns in adata.uns["nhood_adata"].var)
    - model_contrasts: A string vector that defines the contrasts used to perform DA testing
    - subset_samples: subset of samples (columns in `adata.uns["nhood_adata"].X`) to use for the test
    - add_intercept: whether to include an intercept in the model. If False, this is equivalent to adding + 0 in the design formula.
    When model_contrasts is specified, this is set to False by default. 
    '''
    # Set up rpy2 to run edgeR
    rpy2.robjects.numpy2ri.activate()
    rpy2.robjects.pandas2ri.activate()
    edgeR = _try_import_bioc_library("edgeR")
    limma = _try_import_bioc_library("limma")
    stats = importr("stats")
    base = importr("base")

    nhood_adata = adata.uns["nhood_adata"]
    covariates = [x.strip(" ") for x in set(
        re.split('\\+|\\*', design.lstrip("~ ")))]

    # Add covariates used for testing to nhood_adata.var
    sample_col = nhood_adata.uns["sample_col"]
    try:
        nhoods_var = adata.obs[covariates + [sample_col]].drop_duplicates()
    except KeyError:
        missing_cov = [
            x for x in covariates if x not in nhood_adata.var.columns]
        raise KeyError(
            'Covariates {c} are not columns in adata.obs'.format(
                c=" ".join(missing_cov))
        )
    # N.B. This might need some type adjustment!!
    nhoods_var = nhoods_var[covariates + [sample_col]]
    nhoods_var.index = nhoods_var[sample_col].astype("str")

    try:
        assert nhoods_var.loc[nhood_adata.var_names].shape[0] == len(
            nhood_adata.var_names)
    except:
        raise ValueError(
            "Covariates cannot be unambiguously assigned to each sample -- each sample value should match a single covariate value")
    nhood_adata.var = nhoods_var.loc[nhood_adata.var_names]
    # Get design dataframe
    try:
        design_df = nhood_adata.var[covariates]
    except KeyError:
        missing_cov = [
            x for x in covariates if x not in nhood_adata.var.columns]
        raise KeyError(
            'Covariates {c} are not columns in adata.uns["nhood_adata"].var'.format(
                c=" ".join(missing_cov))
        )

    # Get count matrix
    count_mat = nhood_adata.X.toarray()
    lib_size = count_mat.sum(0)

    # Filter out samples with zero counts
    keep_smp = lib_size > 0

    # Subset samples
    if subset_samples is not None:
        keep_smp = keep_smp & nhood_adata.var_names.isin(subset_samples)
    
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
    res.index = nhood_adata.obs_names[keep_nhoods]
    if any([x in nhood_adata.obs.columns for x in res.columns]):
        nhood_adata.obs = nhood_adata.obs.drop(res.columns, axis=1)
    nhood_adata.obs = pd.concat([nhood_adata.obs, res], axis=1)

    # Run Graph spatial FDR correction
    _graph_spatialFDR(adata)


def _graph_spatialFDR(adata):
    '''
    FDR correction weighted on inverse of connectivity of neighbourhoods.
    The distance to the k-th nearest neighbor is used as a measure of connectivity.
    '''

    # use 1/connectivity as the weighting for the weighted BH adjustment from Cydar
    w = 1/adata.uns["nhood_adata"].obs['kth_distance']
    w[np.isinf(w)] = 0

    # Computing a density-weighted q-value.
    pvalues = adata.uns["nhood_adata"].obs["PValue"]
    keep_nhoods = ~pvalues.isna()  #  Filtering in case of test on subset of nhoods
    o = pvalues[keep_nhoods].argsort()
    pvalues = pvalues[keep_nhoods][o]
    w = w[keep_nhoods][o]

    adjp = pd.Series(np.zeros(shape=pvalues.shape),
                     index=pvalues.index)
    adjp = (sum(w)*pvalues/np.cumsum(w))[::-1].cummin()[::-1]
    adjp[adjp > 1] = 1

    ## Store in anndata
    adata.uns["nhood_adata"].obs["SpatialFDR"] = np.nan
    adata.uns["nhood_adata"].obs.loc[keep_nhoods, "SpatialFDR"] = adjp

## -- UTILS -- ##


def _try_import_bioc_library(name):
    try:
        _r_lib = importr(name)
        return _r_lib
    except PackageNotInstalledError:
        raise RuntimeError(
            f"Install Bioconductor library `{name!r}` first as `BiocManager::install({name!r}).`"
        )
