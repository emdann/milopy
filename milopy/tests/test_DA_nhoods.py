import pytest
import scanpy as sc
import numpy as np
import pandas as pd
from milopy.core import make_nhoods
from milopy.core import count_nhoods
from milopy.core import DA_nhoods


@pytest.fixture
def anndata(seed=42):
    adata = sc.datasets.pbmc68k_reduced()
    make_nhoods(adata)

    ## Simulate experimental condition ##
    np.random.seed(seed)
    adata.obs["condition"] = np.random.choice(["ConditionA", "ConditionB"], size=adata.n_obs, p=[0.5,0.5])
    # we simulate differential abundance in NK cells
    DA_cells = adata.obs["louvain"] == "1"
    adata.obs.loc[DA_cells, "condition"] = np.random.choice(["ConditionA", "ConditionB"], size=sum(DA_cells), p=[0.2,0.8])

    ## Simulate replicates ##
    adata.obs["replicate"] = np.random.choice(["R1", "R2", "R3"], size=adata.n_obs)
    adata.obs["sample"] = adata.obs["replicate"] + adata.obs["condition"]
    count_nhoods(adata, sample_col="sample")
    return adata

def test_missing_covariate(anndata):
    adata = anndata.copy() 
    with pytest.raises(KeyError):
        DA_nhoods(adata, design="~ciaone")

def test_non_unique_covariate(anndata):
    adata = anndata.copy() 
    with pytest.raises(ValueError):
        DA_nhoods(adata, design="~phase")
        
## Check that results make sense
def test_pvalues(anndata):
    adata = anndata.copy() 
    DA_nhoods(adata, design="~condition")
    nhood_adata = adata.uns["nhood_adata"]
    min_p, max_p = nhood_adata.obs["PValue"].min(),nhood_adata.obs["PValue"].min()
    assert (min_p >= 0) & (max_p <= 1), "P-values are not between 0 and 1"
    
def test_fdr(anndata):
    adata = anndata.copy() 
    DA_nhoods(adata, design="~condition")
    nhood_adata = adata.uns["nhood_adata"]
    assert np.all(np.round(nhood_adata.obs["PValue"], 10) <= np.round(nhood_adata.obs["SpatialFDR"], 10)), "FDR is higher than uncorrected P-values"