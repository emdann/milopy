import pytest
import scanpy as sc
import numpy as np
from milopy.core import make_nhoods, count_nhoods
from milopy.utils import annotate_nhoods, annotate_nhoods_continuous


@pytest.fixture
def milo_mdata(seed=42):
    milo_mdata = prep_nhood_matrix(seed)
    return milo_mdata


def prep_nhood_matrix(seed):
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
    milo_mdata = count_nhoods(adata, sample_col='sample')
    return milo_mdata

# --- Annotate continuous tests ---

# Test that mean values are within the expected range


def test_nhood_mean_range(milo_mdata):
    annotate_nhoods_continuous(milo_mdata, anno_col='S_score')
    assert milo_mdata['samples'].var['nhood_S_score'].max(
    ) < milo_mdata['cells'].obs['S_score'].max()
    assert milo_mdata['samples'].var['nhood_S_score'].min(
    ) > milo_mdata['cells'].obs['S_score'].min()

# Test that value corresponds to mean


def test_correct_mean(milo_mdata):
    annotate_nhoods_continuous(milo_mdata, anno_col='S_score')
    i = np.random.choice(np.arange(milo_mdata['samples'].n_obs))
    mean_val_nhood = milo_mdata['cells'].obs[milo_mdata['cells'].obsm['nhoods']
                                             [:, i].toarray() == 1]['S_score'].mean()
    assert milo_mdata['samples'].var['nhood_S_score'][i] == pytest.approx(
        mean_val_nhood, 0.0001)


# --- Annotate labels tests ---

# Test that label fractions are in the correct range

def test_nhood_annotation_frac_range(milo_mdata):
    annotate_nhoods(milo_mdata, anno_col='louvain')
    assert milo_mdata['samples'].var['nhood_annotation_frac'].max() <= 1.0
    assert milo_mdata['samples'].var['nhood_annotation_frac'].min() >= 0.0

# Test continuous covariate gives error


def test_nhood_annotation_cont_gives_error(milo_mdata):
    with pytest.raises(ValueError):
        annotate_nhoods(milo_mdata, anno_col='S_score')
