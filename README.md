# milopy
Basic python implementation of Milo for differential abundance testing on KNN graphs, to ease interoperability with [`scanpy`](https://scanpy.readthedocs.io/en/stable/index.html) pipelines for single-cell analysis. See our [preprint](https://www.biorxiv.org/content/10.1101/2020.11.23.393769v1) for details on the statistical framework.

**This implementation is experimental**, for a robust and fully functional implementation see our R package [`miloR`](https://github.com/MarioniLab/miloR).

## Installation

```
git clone https://github.com/emdann/milopy.git
cd milopy
pip install .
```

## Example usage
```
import scanpy as sc
import milopy
import milopy.core as milo
```

Load example dataset
```
adata = sc.datasets.pbmc3k_processed()
```

Simulate experimental condition and replicates
```
## Simulate experimental condition ##
adata.obs["condition"] = np.random.choice(["ConditionA", "ConditionB"], size=adata.n_obs, p=[0.5,0.5])
# we simulate differential abundance in NK cells
DA_cells = adata.obs["louvain"] == "NK cells"
adata.obs.loc[DA_cells, "condition"] = np.random.choice(["ConditionA", "ConditionB"], size=sum(DA_cells), p=[0.2,0.8])

## Simulate replicates ##
adata.obs["replicate"] = np.random.choice(["R1", "R2", "R3"], size=adata.n_obs)
adata.obs["sample"] = adata.obs["replicate"] + adata.obs["condition"]

sc.pl.umap(adata, color=["louvain","condition", "sample"])
```

Test for differential abundance with Milo
```
## Build KNN graph
sc.pp.neighbors(adata, n_neighbors=10)

## Assign cells to neighbourhoods
milo.make_nhoods(adata)

## Count cells from each sample in each nhood
milo.count_nhoods(adata, sample_col="sample")

## Test for differential abundance between conditions
milo.test_nhoods(adata, design="~ condition")

## Check results
milo_results = adata.uns["nhood_adata"].obs
milo_results
```

Visualize results on UMAP embedding
```
milopy.utils.build_nhood_graph(adata)
milopy.plot.plot_nhood_graph(adata, alpha=0.2, min_size=5)
```