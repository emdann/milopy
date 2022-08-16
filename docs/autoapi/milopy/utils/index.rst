:py:mod:`milopy.utils`
======================

.. py:module:: milopy.utils


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   milopy.utils.add_nhood_expression
   milopy.utils.build_nhood_graph
   milopy.utils.add_covariate_to_nhoods_var
   milopy.utils.annotate_nhoods
   milopy.utils.annotate_nhoods_continuous
   milopy.utils.write_milo_adata
   milopy.utils.read_milo_adata



.. py:function:: add_nhood_expression(adata: anndata.AnnData, layer: str = None)

   Calculates the mean expression in neighbourhoods of each feature in `adata.X` or
   `adata.layers[layer]` (if layer is not None).

   Params:
   -------
   - adata: AnnData object
   - layer: which data layer to use as expression matrix (default: None, uses `adata.X`)

   Returns:
   -------
   Updates adata in place to store the matrix of average expression in each neighbourhood in `adata.uns["nhood_adata"].obsm['expr']`


.. py:function:: build_nhood_graph(adata: anndata.AnnData, basis: str = 'X_umap')

   Build graph of neighbourhoods used for visualization of DA results

   Params:
   -------
   - adata: AnnData object
   - basis: string indicating the name of the obsm basis to use to use for layout of neighbourhoods (key in `adata.obsm`)


.. py:function:: add_covariate_to_nhoods_var(adata: anndata.AnnData, new_covariates: List[str])

   Add covariate from adata.obs to adata.uns["nhood_adata"].var


.. py:function:: annotate_nhoods(adata: anndata.AnnData, anno_col: str)

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


.. py:function:: annotate_nhoods_continuous(adata: anndata.AnnData, anno_col: str)

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


.. py:function:: write_milo_adata(adata: anndata.AnnData, filepath: str, **kwargs)

   Save anndata objects after Milo analysis

   Params:
   -----
   - adata: AnnData object with adata.uns["nhood_adata"]
   - filepath: path to h5ad file to save
   - **kwargs: arguments passed to scanpy.write_h5ad 

   Returns:
   -------
   None, saves 2 AnnData objects in h5ad format. The cell x gene AnnData is saved in filepath.
   The nhood x sample AnnData is saved in a separate object (location is stored in adata.uns['nhood_adata_filepath'])


.. py:function:: read_milo_adata(filepath: str, **kwargs) -> anndata.AnnData

   Read AnnData objects stored after Milo analysis

   Params:
   ------
   - filepath: path to h5ad file storing cell x gene AnnData object
   - **kwargs: additional arguments passed to scanpy.read_h5ad

   Returns:
   -------
   - AnnData object storing milo slots (adata.obsm['nhoods'], adata.uns['nhood_adata'])


