:py:mod:`milopy.plot`
=====================

.. py:module:: milopy.plot


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   milopy.plot.plot_nhood_graph
   milopy.plot.plot_nhood
   milopy.plot.plot_DA_beeswarm
   milopy.plot._get_palette_adata
   milopy.plot.plot_nhood_counts_by_cond



.. py:function:: plot_nhood_graph(adata: anndata.AnnData, alpha: float = 0.1, min_logFC: float = 0, min_size: int = 10, plot_edges: bool = False, title: str = 'DA log-Fold Change', **kwargs)

   Visualize DA results on abstracted graph (wrapper around sc.pl.embedding)

   Params:
   -------
   - adata: AnnData object
   - alpha: significance threshold
   - min_logFC: minimum absolute log-Fold Change to show results (default: 0, show all significant neighbourhoods)
   - min_size: minimum size of nodes in visualization (default: 10)
   - plot_edges: boolean indicating if edges for neighbourhood overlaps whould be plotted (default: False)
   - title: plot title (default: 'DA log-Fold Change')
   - **kwargs: other arguments to pass to scanpy.pl.embedding


.. py:function:: plot_nhood(adata, ix, basis='X_umap')

   Visualize cells in a neighbourhood


.. py:function:: plot_DA_beeswarm(adata: anndata.AnnData, anno_col: str = 'nhood_annotation', alpha: float = 0.1, subset_nhoods: List = None)

   Plot beeswarm plot of logFC against nhood labels

   Params:
   -------
   - adata: AnnData object
   - anno_col: column in adata.uns['nhood_adata'].obs to use as annotation
   - alpha: significance threshold
   - subset_nhoods: list of nhoods to plot (default: None, plot all nhoods)


.. py:function:: _get_palette_adata(adata, obs_col)


.. py:function:: plot_nhood_counts_by_cond(adata: anndata.AnnData, test_var: str, subset_nhoods: List = None, log_counts: bool = False)

   Plot boxplot of cell numbers vs condition of interest

   Params:
   ------
   - adata: anndata object storing neighbourhood information in adata.uns
   - test_var: string, name of column in adata.obs storing condition of interest (y-axis for boxplot)
   - subset_nhoods: list of obs_names for neighbourhoods to include in plot (default: None, plot all nhoods)
   - log_counts: boolean, whether to plot log1p of cell counts (default: False)


