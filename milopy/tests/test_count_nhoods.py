import pytest
import scanpy as sc
import numpy as np
import pandas as pd
from milopy.core import make_nhoods
from milopy.core import count_nhoods
       
adata = sc.datasets.pbmc68k_reduced()
make_nhoods(adata)
    
## Test that
def test_sample_values():
    ## Extract cells of one nhood
    nh=1
    sample_col="phase"
    count_nhoods(adata, sample_col=sample_col)
    nh_cells = adata.obsm["nhoods"][:,nh].nonzero()[0]

    ## Value count the sample composition
    top_a = adata.obs.iloc[nh_cells].value_counts(sample_col).values.ravel()

    ## Check it matches the one calculated
    df = pd.DataFrame(adata.uns["nhood_adata"].X[nh,:].toarray()).T
    df.index = adata.uns["nhood_adata"].var_names
    top_b = df.sort_values(0, ascending=False).values.ravel()
    assert all((top_b - top_a) == 0), 'The counts for samples in adata.uns["nhood_adata"] does not match'
    
def test_sample_order():    
    ## Extract cells of one nhood
    nh=1
    sample_col="phase"
    count_nhoods(adata, sample_col=sample_col)
    nh_cells = adata.obsm["nhoods"][:,nh].nonzero()[0]

    ## Value count the sample composition
    top_a = adata.obs.iloc[nh_cells].value_counts(sample_col).index[0]

    ## Check it matches the one calculated
    df = pd.DataFrame(adata.uns["nhood_adata"].X[nh,:].toarray()).T
    df.index = adata.uns["nhood_adata"].var_names
    top_b = df.sort_values(0, ascending=False).index[0]

    assert top_a==top_b, 'The order of samples in adata.uns["nhood_adata"] does not match'
 