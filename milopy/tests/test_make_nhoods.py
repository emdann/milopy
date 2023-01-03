### TESTS ###
import pytest
import scanpy as sc
import numpy as np
from milopy.core import make_nhoods

adata_example = sc.datasets.pbmc68k_reduced()

# Test there are less nhoods than n_obs*prop


def test_number_of_nhoods():
    p = 0.1
    make_nhoods(adata_example, prop=p)
    assert adata_example.obsm["nhoods"].shape[1] <= int(
        np.round(adata_example.n_obs * p))

# Test that the smallest neighbourhood is bigger or equal than
# the smallest number of neighbors in the KNN graph


def test_nhood_sizes():
    make_nhoods(adata_example)
    knn_graph = adata_example.obsp["connectivities"]
    knn_graph[knn_graph != 0] = 1
    assert knn_graph.sum(0).min() <= adata_example.obsm["nhoods"].sum(0).min()

## Test that the right KNN graph is used if specified ##

def test_neighbors_key():
    k = adata_example.uns['neighbors']['params']['n_neighbors']
    test_k = 5
    assert test_k < k
    sc.pp.neighbors(adata_example, n_neighbors=test_k, key_added='test')
    make_nhoods(adata_example, neighbors_key='test')
    smallest_size = adata_example.obsm['nhoods'].toarray().sum(0).min()
    assert smallest_size < k
    assert adata_example.uns["nhood_neighbors_key"] == 'test'
