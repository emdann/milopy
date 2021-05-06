import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple
from anndata import AnnData

import matplotlib.pyplot as plt
import seaborn as sns

def plot_nhood_graph(
    adata, 
    alpha=0.1,
    min_size=10, 
    **kwargs
    ):
    '''
    Visualize DA results on abstracted graph (wrapper around sc.pl.embedding)
    
    - adata: AnnData object
    - alpha: significance threshold
    - min_size: minimum size of nodes in visualization (default: 10)
    - **kwargs: other arguments to pass to scanpy.pl.embedding
    '''
    nhood_adata = adata.uns["nhood_adata"]
    
    vmax = np.max([nhood_adata.obs["logFC"].max(), abs(nhood_adata.obs["logFC"].min())])
    vmin = - vmax

    nhood_adata.obs["graph_color"] = nhood_adata.obs["logFC"]
    nhood_adata.obs.loc[nhood_adata.obs["SpatialFDR"] > alpha, "graph_color"] = np.nan
    nhood_adata.obs["graph_color"].max()

    sc.pl.embedding(nhood_adata, "X_milo_graph", 
                    color="graph_color", cmap="RdBu",
                    size=adata.uns["nhood_adata"].obs["Nhood_size"]*min_size, 
                    edges=True, neighbors_key="nhood",
                    # edge_width = 
                    frameon=False,
                    vmax=vmax, vmin=vmin,
                    title="DA log-Fold Change"
                   )
    
def plot_nhood(adata, ix, basis="X_umap"):    
    '''
    Visualize cells in a neighbourhood
    '''
    adata.obs["Nhood"] = adata.obsm["nhoods"][:,ix].toarray().ravel()
    sc.pl.embedding(adata, basis, color="Nhood", size=30, title="Nhood" + str(ix))