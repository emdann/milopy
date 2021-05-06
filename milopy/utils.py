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

    ## Get gene expression matrix
    if layer is None:
        X = adata.X.value
        expr_id = "expr" 
    else:
        X = adata.layers[layer]
        expr_id = "expr_" + layer

    ## Aggregate over nhoods -- taking the mean
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
    ## Add embedding positions
    adata.uns["nhood_adata"].obsm["X_milo_graph"] = adata[adata.obs["nhood_ixs_refined"]==1].obsm[use_rep]
    ## Add nhood size
    adata.uns["nhood_adata"].obs["Nhood_size"] = np.array(adata.obsm["nhoods"].sum(0)).flatten()
    ## Add adjacency graph
    adata.uns["nhood_adata"].obsp["nhood_connectivities"] = adata.obsm["nhoods"].T.dot(adata.obsm["nhoods"])
    adata.uns["nhood_adata"].uns["nhood"] = {"connectivities_key":"nhood_connectivities", "distances_key":""}

