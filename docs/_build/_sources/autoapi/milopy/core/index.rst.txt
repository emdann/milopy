:py:mod:`milopy.core`
=====================

.. py:module:: milopy.core


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   milopy.core.make_nhoods
   milopy.core.count_nhoods
   milopy.core.DA_nhoods
   milopy.core._graph_spatialFDR
   milopy.core._try_import_bioc_library



.. py:function:: make_nhoods(adata: anndata.AnnData, neighbors_key: str = None, prop=0.1, seed: int = 42)

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


.. py:function:: count_nhoods(adata: anndata.AnnData, sample_col: str)

   - adata
   - sample_col: string, column in adata.obs that contains sample information 
   (what should be in the columns of the nhoodCount matrix)

   Returns: None
   Updated adata.uns slot to contain adata.uns["nhood_adata"], where:
   - adata.uns["nhood_adata"].obs_names are neighbourhoods
   - adata.uns["nhood_adata"].var_names are samples
   - adata.uns["nhood_adata"].X is the matrix counting the number of cells from each
   sample in each neighbourhood


.. py:function:: DA_nhoods(adata, design, model_contrasts=None, subset_samples=None, add_intercept=True)

   This will perform differential neighbourhood abundance testing (using edgeR under the hood)
   - adata
   - design: formula (terms should be columns in adata.uns["nhood_adata"].var)
   - model_contrasts: A string vector that defines the contrasts used to perform DA testing
   - subset_samples: subset of samples (columns in `adata.uns["nhood_adata"].X`) to use for the test
   - add_intercept: whether to include an intercept in the model. If False, this is equivalent to adding + 0 in the design formula.
   When model_contrasts is specified, this is set to False by default. 


.. py:function:: _graph_spatialFDR(adata, neighbors_key=None)

   FDR correction weighted on inverse of connectivity of neighbourhoods.
   The distance to the k-th nearest neighbor is used as a measure of connectivity.


.. py:function:: _try_import_bioc_library(name)


