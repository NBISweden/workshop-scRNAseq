---
description: Reduce high-dimensional gene expression data from
subtitle:  SCANPY TOOLKIT
title:  Dimensionality Reduction
---

<div>

> **Note**
>
> Code chunks run Python commands unless it starts with `%%bash`, in
> which case, those chunks run shell commands.

</div>

## Data preparation

First, let's load all necessary libraries and the QC-filtered dataset
from the previous step.

``` {python}
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 3
#sc.logging.print_versions()

sc.settings.set_figure_params(dpi=80)
adata = sc.read_h5ad('data/results/scanpy_qc_filtered_covid.h5ad')
adata
```

Before variable gene selection we need to normalize and log transform
the data. Then store the full matrix in the `raw` slot before doing
variable gene selection.

``` {python}
# normalize to depth 10 000
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

# log transform
sc.pp.log1p(adata)

# store normalized counts in the raw slot, 
# we will subset adata.X for variable genes, but want to keep all genes matrix as well.
adata.raw = adata

adata
```

## Feature selection

Next, we first need to define which features/genes are important in our
dataset to distinguish cell types. For this purpose, we need to find
genes that are highly variable across cells, which in turn will also
provide a good separation of the cell clusters.

``` {python}
# compute variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print("Highly variable genes: %d"%sum(adata.var.highly_variable))

#plot variable genes
sc.pl.highly_variable_genes(adata)

# subset for variable genes in the dataset
adata = adata[:, adata.var['highly_variable']]
```

## Z-score transformation

Now that the data is prepared, we now proceed with PCA. Since each gene
has a different expression level, it means that genes with higher
expression values will naturally have higher variation that will be
captured by PCA. This means that we need to somehow give each gene a
similar weight when performing PCA (see below). The common practice is
to center and scale each gene before performing PCA. This exact scaling
is called Z-score normalization it is very useful for PCA, clustering
and plotting heatmaps. Additionally, we can use regression to remove any
unwanted sources of variation from the dataset, such as `cell cycle`,
`sequencing depth`, `percent mitochondria`. This is achieved by doing a
generalized linear regression using these parameters as co-variates in
the model. Then the residuals of the model are taken as the *regressed
data*. Although perhaps not in the best way, batch effect regression can
also be done here. By default variables are scaled in the PCA step and
is not done separately. But it could be achieved by running the commands
below:

``` {python}
#run this line if you get the "AttributeError: swapaxes not found" 
# adata = adata.copy()

# regress out unwanted variables
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# scale data, clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)
```

## PCA

Performing PCA has many useful applications and interpretations, which
much depends on the data used. In the case of life sciences, we want to
segregate samples based on gene expression patterns in the data.

To run PCA, you can use the function `pca()`.

``` {python}
sc.tl.pca(adata, svd_solver='arpack')
```

We then plot the first principal components.

``` {python}
# plot more PCS
sc.pl.pca(adata, color='sample', components = ['1,2','3,4','5,6','7,8'], ncols=2)
```

To identify genes that contribute most to each PC, one can retrieve the
loading matrix information.

``` {python}
#Plot loadings
sc.pl.pca_loadings(adata, components=[1,2,3,4,5,6,7,8])

# OBS! only plots the positive axes genes from each PC!!
```

The function to plot loading genes only plots genes on the positive
axes. Instead plot as a heatmaps, with genes on both positive and
negative side, one per pc, and plot their expression amongst cells
ordered by their position along the pc.

``` {python}
# adata.obsm["X_pca"] is the embeddings
# adata.uns["pca"] is pc variance
# adata.varm['PCs'] is the loadings

genes = adata.var['gene_ids']

for pc in [1,2,3,4]:
    g = adata.varm['PCs'][:,pc-1]
    o = np.argsort(g)
    sel = np.concatenate((o[:10],o[-10:])).tolist()
    emb = adata.obsm['X_pca'][:,pc-1]
    # order by position on that pc
    tempdata = adata[np.argsort(emb),]
    sc.pl.heatmap(tempdata, var_names = genes[sel].index.tolist(), groupby='predicted_doublets', swap_axes = True, use_raw=False)
```

We can also plot the amount of variance explained by each PC.

``` {python}
sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 50)
```

Based on this plot, we can see that the top 8 PCs retain a lot of
information, while other PCs contain progressively less. However, it is
still advisable to use more PCs since they might contain information
about rare cell types (such as platelets and DCs in this dataset)

## tSNE

We can now run [BH-tSNE](https://arxiv.org/abs/1301.3342).

``` {python}
sc.tl.tsne(adata, n_pcs = 30)
```

We can now plot the tSNE colored per dataset. We can clearly see the
effect of batches present in the dataset.

``` {python}
sc.pl.tsne(adata, color='sample')
```

## UMAP

The UMAP implementation in SCANPY uses a neighborhood graph as the
distance matrix, so we need to first calculate the graph.

``` {python}
sc.pp.neighbors(adata, n_pcs = 30, n_neighbors = 20)
```

We can now run [UMAP](https://arxiv.org/abs/1802.03426) for cell
embeddings.

``` {python}
sc.tl.umap(adata)
sc.pl.umap(adata, color='sample')
```

UMAP is plotted colored per dataset. Although less distinct as in the
tSNE, we still see quite an effect of the different batches in the data.
UMAP is not limited by the number of dimensions the data can be reduced
into (unlike tSNE). We can simply reduce the dimensions altering the
`n.components` parameter.

``` {python}
#run with 10 components, save to a new object so that the umap with 2D is not overwritten.
umap10 = sc.tl.umap(adata, n_components=10, copy=True)
fig, axs = plt.subplots(1, 3, figsize=(10,4),constrained_layout=True)

sc.pl.umap(adata, color='sample',  title="UMAP", show=False, ax=axs[0])
sc.pl.umap(umap10, color='sample', title="UMAP10", show=False, ax=axs[1], components=['1,2'])
sc.pl.umap(umap10, color='sample', title="UMAP10", show=False, ax=axs[2], components=['3,4'])

# we can also plot the umap with neighbor edges
sc.pl.umap(adata, color='sample', title="UMAP", edges=True)
```

We can now plot PCA, UMAP and tSNE side by side for comparison. Here, we
can conclude that our dataset contains a batch effect that needs to be
corrected before proceeding to clustering and differential gene
expression analysis.

TODO: \[pca, tsne and umap plots side by side\]

<div>

> **Discuss**
>
> We have now done Variable gene selection, PCA and UMAP with the
> settings we chose. Test a few different ways of selecting variable
> genes, number of PCs for UMAP and check how it influences your
> embedding.

</div>

## Genes of interest

Let's plot some marker genes for different cell types onto the
embedding.

  Markers                    Cell Type
  -------------------------- -------------------
  CD3E                       T cells
  CD3E CD4                   CD4+ T cells
  CD3E CD8A                  CD8+ T cells
  GNLY, NKG7                 NK cells
  MS4A1                      B cells
  CD14, LYZ, CST3, MS4A7     CD14+ Monocytes
  FCGR3A, LYZ, CST3, MS4A7   FCGR3A+ Monocytes
  FCER1A, CST3               DCs

``` {python}
sc.pl.umap(adata, color=["CD3E", "CD4", "CD8A", "GNLY","NKG7", "MS4A1","CD14","LYZ","CST3","MS4A7","FCGR3A"])
```

The default is to plot gene expression in the normalized and
log-transformed data. You can also plot it on the scaled and corrected
data by using `use_raw=False`. However, not all of these genes are
included in the variable gene set so we first need to filter them.

``` {python}
genes  = ["CD3E", "CD4", "CD8A", "GNLY","NKG7", "MS4A1","CD14","LYZ","CST3","MS4A7","FCGR3A"]
var_genes = adata.var.highly_variable
var_genes.index[var_genes]
varg = [x for x in genes if x in var_genes.index[var_genes]]
sc.pl.umap(adata, color=varg, use_raw=False)
```

<div>

> **Discuss**
>
> Select some of your dimensionality reductions and plot some of the QC
> stats that were calculated in the previous lab. Can you see if some of
> the separation in your data is driven by quality of the cells?

</div>

## Save data

We can finally save the object for use in future steps.

``` {python}
import os
path = "./data/results/scanpy_dr_covid.h5ad"
if not os.path.exists(path):
    adata.write_h5ad(path)
```

``` {python}
print(adata.X.shape)
print(adata.raw.X.shape)
```

## Session info

``` {python}
sc.logging.print_versions()
```
