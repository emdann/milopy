import pytest
import scanpy as sc
import numpy as np
import pandas as pd
from milopy.core import make_nhoods
from milopy.core import count_nhoods
from milopy.core import DA_nhoods


@pytest.fixture
def milo_mdata(seed=42):
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
    milo_mdata = count_nhoods(adata, sample_col="sample")
    return milo_mdata


def test_missing_covariate(milo_mdata):
    milo_mdata = milo_mdata.copy()
    with pytest.raises(KeyError):
        DA_nhoods(milo_mdata, design="~ciaone")


def test_non_unique_covariate(milo_mdata):
    milo_mdata = milo_mdata.copy()
    with pytest.raises(ValueError):
        DA_nhoods(milo_mdata, design="~phase")

# Check that results make sense


def test_pvalues(milo_mdata):
    milo_mdata = milo_mdata.copy()
    DA_nhoods(milo_mdata, design="~condition")
    sample_adata = milo_mdata['samples'].copy()
    min_p, max_p = sample_adata.var["PValue"].min(
    ), sample_adata.var["PValue"].max()
    assert (min_p >= 0) & (max_p <= 1), "P-values are not between 0 and 1"


def test_fdr(milo_mdata):
    milo_mdata = milo_mdata.copy()
    DA_nhoods(milo_mdata, design="~condition")
    sample_adata = milo_mdata['samples'].copy()
    assert np.all(np.round(sample_adata.var["PValue"], 10) <= np.round(
        sample_adata.var["SpatialFDR"], 10)), "FDR is higher than uncorrected P-values"

# Test specifying contrasts


def test_default_contrast(milo_mdata):
    milo_mdata = milo_mdata.copy()
    adata = milo_mdata['cells'].copy()
    adata.obs['condition'] = adata.obs['condition'].astype(
        'category').cat.reorder_categories(['ConditionA', 'ConditionB'])
    DA_nhoods(milo_mdata, design='~condition')
    default_results = milo_mdata['samples'].var.copy()
    DA_nhoods(milo_mdata, design='~condition',
              model_contrasts='conditionConditionB-conditionConditionA')
    contr_results = milo_mdata['samples'].var.copy()

    assert np.corrcoef(contr_results['SpatialFDR'], default_results['SpatialFDR'])[
        0, 1] > 0.99
    assert np.corrcoef(contr_results['logFC'],
                       default_results['logFC'])[0, 1] > 0.99
