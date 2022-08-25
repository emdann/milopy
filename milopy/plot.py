import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from typing import Union, Optional, Sequence, Any, Mapping, List, Tuple
from anndata import AnnData
from mudata import MuData

import matplotlib.pyplot as plt
import seaborn as sns


def plot_nhood_graph(
    milo_mdata: MuData,
    alpha: float = 0.1,
    min_logFC: float = 0,
    min_size: int = 10,
    plot_edges: bool = False,
    title: str = "DA log-Fold Change",
    **kwargs
):
    '''
    Visualize DA results on abstracted graph (wrapper around sc.pl.embedding)

    Params:
    -------
    - milo_mdata: MuData object
    - alpha: significance threshold
    - min_logFC: minimum absolute log-Fold Change to show results (default: 0, show all significant neighbourhoods)
    - min_size: minimum size of nodes in visualization (default: 10)
    - plot_edges: boolean indicating if edges for neighbourhood overlaps whould be plotted (default: False)
    - title: plot title (default: 'DA log-Fold Change')
    - **kwargs: other arguments to pass to scanpy.pl.embedding
    '''
    nhood_adata = milo_mdata['samples'].T.copy()

    if "Nhood_size" not in nhood_adata.obs.columns:
        raise KeyError(
            'Cannot find "Nhood_size" column in adata.uns["nhood_adata"].obs -- \
                please run milopy.utils.build_nhood_graph(adata)'
        )

    nhood_adata.obs["graph_color"] = nhood_adata.obs["logFC"]
    nhood_adata.obs.loc[nhood_adata.obs["SpatialFDR"]
                        > alpha, "graph_color"] = np.nan
    nhood_adata.obs["abs_logFC"] = abs(nhood_adata.obs["logFC"])
    nhood_adata.obs.loc[nhood_adata.obs["abs_logFC"]
                        < min_logFC, "graph_color"] = np.nan

    # Plotting order - extreme logFC on top
    nhood_adata.obs.loc[nhood_adata.obs["graph_color"].isna(),
                        "abs_logFC"] = np.nan
    ordered = nhood_adata.obs.sort_values(
        'abs_logFC', na_position='first').index
    nhood_adata = nhood_adata[ordered]

    vmax = np.max([nhood_adata.obs["graph_color"].max(),
                  abs(nhood_adata.obs["graph_color"].min())])
    vmin = - vmax

    sc.pl.embedding(nhood_adata, "X_milo_graph",
                    color="graph_color", cmap="RdBu_r",
                    size=nhood_adata.obs["Nhood_size"]*min_size,
                    edges=plot_edges, neighbors_key="nhood",
                    sort_order=False,
                    frameon=False,
                    vmax=vmax, vmin=vmin,
                    title=title,
                    **kwargs
                    )


def plot_nhood(adata, ix, basis="X_umap"):
    '''
    Visualize cells in a neighbourhood

    Params:
    ------
    - adata: AnnData object, storing neighbourhood assignments in `adata.obsm['nhoods']`
    - ix: index of neighbourhood to visualize
    - basis: embedding to use for visualization (default: 'X_umap')
    '''
    adata.obs["Nhood"] = adata.obsm["nhoods"][:, ix].toarray().ravel()
    sc.pl.embedding(adata, basis, color="Nhood",
                    size=30, title="Nhood" + str(ix))

# Plot nhood beeswarm


def plot_DA_beeswarm(
    milo_mdata: MuData,
    anno_col: str = 'nhood_annotation',
    alpha: float = 0.1,
    subset_nhoods: List = None
):
    '''
    Plot beeswarm plot of logFC against nhood labels

    Params:
    -------
    - milo_mdata: MuData object
    - anno_col: column in adata.uns['nhood_adata'].obs to use as annotation
    - alpha: significance threshold
    - subset_nhoods: list of nhoods to plot (default: None, plot all nhoods)
    '''
    try:
        nhood_adata = milo_mdata['samples'].T.copy()
    except KeyError:
        raise KeyError(
            "milo_mdata should be a MuData object with two slots: 'cells' and 'samples' - please run milopy.count_nhoods(adata) first"
        )

    if subset_nhoods is not None:
        nhood_adata = nhood_adata[subset_nhoods]

    try:
        nhood_adata.obs[anno_col]
    except KeyError:
        raise KeyError(
            'Cannot find {a} in adata.uns["nhood_adata"].obs -- \
                please run milopy.utils.annotate_nhoods(adata, anno_col) first'.format(a=anno_col)
        )

    try:
        nhood_adata.obs["logFC"]
    except KeyError:
        raise KeyError(
            'Cannot find `logFC` in adata.uns["nhood_adata"].obs -- \
                please run milopy.core.DA_nhoods(adata) first'
        )

    sorted_annos = nhood_adata.obs[[anno_col, "logFC"]].\
        groupby(anno_col).\
        median().\
        sort_values("logFC", ascending=True).index

    anno_df = nhood_adata.obs[[anno_col, "logFC", "SpatialFDR"]].copy()
    anno_df['is_signif'] = anno_df['SpatialFDR'] < alpha
    # anno_df['across_organs'] = ["Significant across organs" if x else "" for x in (keep_nhoods & signif_nhoods)]
    anno_df = anno_df[anno_df[anno_col] != "nan"]

    try:
        anno_palette = _get_palette_adata(
            milo_mdata['cells'], nhood_adata.uns['annotation_obs'])
        sns.violinplot(data=anno_df, y=anno_col, x="logFC", order=sorted_annos,
                       size=190, inner=None, orient="h",
                       palette=anno_palette,
                       linewidth=0,
                       scale="width")
    except:
        sns.violinplot(data=anno_df, y=anno_col, x="logFC", order=sorted_annos,
                       size=190, inner=None, orient="h",
                       linewidth=0,
                       scale="width")
    sns.stripplot(data=anno_df, y=anno_col, x="logFC", order=sorted_annos, size=2,
                  hue='is_signif', palette=['grey', 'black'],
                  orient="h", alpha=0.5)
    plt.legend(loc='upper left', title=f'< {int(alpha*100)}% SpatialFDR',
               bbox_to_anchor=(1, 1), frameon=False)
    plt.axvline(x=0, ymin=0, ymax=1, color="black", linestyle="--")


def _get_palette_adata(adata, obs_col):
    return(dict(zip(adata.obs[obs_col].cat.categories, adata.uns[f'{obs_col}_colors'])))


### Plot boxplot of cell numbers for QC ###
def plot_nhood_counts_by_cond(
        milo_mdata: MuData,
        test_var: str,
        subset_nhoods: List = None,
        log_counts: bool = False):
    '''
    Plot boxplot of cell numbers vs condition of interest

    Params:
    ------
    - milo_mdata: MuData object storing cell level and nhood level information
    - test_var: string, name of column in adata.obs storing condition of interest (y-axis for boxplot)
    - subset_nhoods: list of obs_names for neighbourhoods to include in plot (default: None, plot all nhoods)
    - log_counts: boolean, whether to plot log1p of cell counts (default: False)
    '''
    try:
        nhood_adata = milo_mdata['samples'].T.copy()
    except KeyError:
        raise KeyError(
            "milo_mdata should be a MuData object with two slots: 'cells' and 'samples' - please run milopy.count_nhoods(adata) first"
        )

    if subset_nhoods is None:
        subset_nhoods = nhood_adata.obs_names

    pl_df = pd.DataFrame(nhood_adata[subset_nhoods].X.A, columns=nhood_adata.var_names).melt(
        var_name=nhood_adata.uns['sample_col'], value_name='n_cells')
    pl_df = pd.merge(pl_df, nhood_adata.var)
    pl_df['log_n_cells'] = np.log1p(pl_df['n_cells'])
    if not log_counts:
        sns.boxplot(data=pl_df, x=test_var, y='n_cells', color='lightblue')
        sns.stripplot(data=pl_df, x=test_var, y='n_cells', color='black', s=3)
        plt.ylabel('# cells')
    else:
        sns.boxplot(data=pl_df, x=test_var, y='log_n_cells', color='lightblue')
        sns.stripplot(data=pl_df, x=test_var,
                      y='log_n_cells', color='black', s=3)
        plt.ylabel('log(# cells + 1)')

    plt.xticks(rotation=90)
    plt.xlabel(test_var)
