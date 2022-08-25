import pytest
import scanpy as sc
import numpy as np
from milopy.core import make_nhoods, count_nhoods
from milopy.utils import add_nhood_expression


@pytest.fixture
def mdata(seed=42):
    adata = sc.datasets.pbmc3k()
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    make_nhoods(adata)

    ## Simulate experimental condition ##
    np.random.seed(seed)
    adata.obs["condition"] = np.random.choice(
        ["ConditionA", "ConditionB"], size=adata.n_obs, p=[0.2, 0.8])
    # we simulate differential abundance in NK cells
    DA_cells = adata.obs["leiden"] == "1"
    adata.obs.loc[DA_cells, "condition"] = np.random.choice(
        ["ConditionA", "ConditionB"], size=sum(DA_cells), p=[0.2, 0.8])

    ## Simulate replicates ##
    adata.obs["replicate"] = np.random.choice(
        ["R1", "R2", "R3"], size=adata.n_obs)
    adata.obs["sample"] = adata.obs["replicate"] + adata.obs["condition"]
    milo_mdata = count_nhoods(adata, sample_col="sample")
    return milo_mdata

## Test that dimensions are correct ##


def test_nhood_mean_range(mdata):
    milo_mdata = mdata.copy()
    add_nhood_expression(milo_mdata)
    assert milo_mdata['samples'].varm['expr'].shape[1] == milo_mdata['cells'].n_vars

## Test the result is the mean ##


def test_nhood_mean_range(mdata):
    milo_mdata = mdata.copy()
    add_nhood_expression(milo_mdata)
    nhood_ix = 10
    nhood_gex = milo_mdata['samples'].varm['expr'][nhood_ix,
                                                   :].toarray().ravel()
    nhood_cells = milo_mdata['cells'].obs_names[milo_mdata['cells'].obsm['nhoods']
                                                [:, nhood_ix].toarray().ravel() == 1]
    mean_gex = np.array(milo_mdata['cells']
                        [nhood_cells].X.mean(axis=0)).ravel()
    nhood_gex == pytest.approx(mean_gex, 0.0001)
