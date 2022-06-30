import pytest
from milopy.utils import write_milo_adata
import scanpy as sc
import numpy as np
import pandas as pd
from milopy.core import make_nhoods
from milopy.core import count_nhoods
from milopy.utils import write_milo_adata


@pytest.fixture
def anndata(seed=42):
    adata = sc.datasets.pbmc68k_reduced()
    make_nhoods(adata)

    ## Simulate experimental condition ##
    np.random.seed(seed)
    adata.obs["condition"] = np.random.choice(
        ["ConditionA", "ConditionB"], size=adata.n_obs, p=[0.5, 0.5])
    # we simulate differential abundance in NK cells
    DA_cells = adata.obs["louvain"] == "1"
    adata.obs.loc[DA_cells, "condition"] = np.random.choice(
        ["ConditionA", "ConditionB"], size=sum(DA_cells), p=[0.2, 0.8])

    ## Simulate replicates ##
    adata.obs["replicate"] = np.random.choice(
        ["R1", "R2", "R3"], size=adata.n_obs)
    adata.obs["sample"] = adata.obs["replicate"] + adata.obs["condition"]
    count_nhoods(adata, sample_col="sample")
    return adata


def test_raise_error(anndata):
    adata = anndata.copy()
    del adata.uns['nhood_adata']
    with pytest.raises(KeyError):
        write_milo_adata(adata, design="~ciaone")
