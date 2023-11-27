---
description: Grouping individual cells with similar gene expression
subtitle:  SEURAT TOOLKIT
title:  Clustering
---

<div>

> **Note**
>
> Code chunks run R commands unless otherwise specified.

</div>

In this tutorial we will continue the analysis of the integrated
dataset. We will use the integrated PCA to perform the clustering. First
we will construct a $k$-nearest neighbor graph in order to perform a
clustering on the graph. We will also show how to perform hierarchical
clustering and k-means clustering on PCA space.

Let's first load all necessary libraries and also the integrated dataset
from the previous step.

``` {r}

suppressPackageStartupMessages({
  library(Seurat)
  library(cowplot)
  library(ggplot2)
  library(pheatmap)
  library(rafalib)
  library(clustree)
})

alldata <- readRDS("data/results/covid_qc_dr_int.rds")
```

## Graph clustering

The procedure of clustering on a Graph can be generalized as 3 main
steps: 1) Build a kNN graph from the data 2) Prune spurious connections
from kNN graph (optional step). This is a SNN graph. 3) Find groups of
cells that maximizes the connections within the group compared other
groups.

### Building kNN / SNN graph

The first step into graph clustering is to construct a k-nn graph, in
case you don't have one. For this, we will use the PCA space. Thus, as
done for dimensionality reduction, we will use ony the top *N* PCA
dimensions for this purpose (the same used for computing UMAP / tSNE).

As we can see above, the **Seurat** function `FindNeighbors()` already
computes both the KNN and SNN graphs, in which we can control the
minimal percentage of shared neighbours to be kept. See `?FindNeighbors`
for additional options.

``` {r}
# check that CCA is still the active assay
alldata@active.assay

alldata <- FindNeighbors(alldata, dims = 1:30, k.param = 60, prune.SNN = 1/15)

# check the names for graphs in the object.
names(alldata@graphs)
```

We can take a look at the kNN graph. It is a matrix where every
connection between cells is represented as $1$s. This is called a
**unweighted** graph (default in Seurat). Some cell connections can
however have more importance than others, in that case the scale of the
graph from $0$ to a maximum distance. Usually, the smaller the distance,
the closer two points are, and stronger is their connection. This is
called a **weighted** graph. Both weighted and unweighted graphs are
suitable for clustering, but clustering on unweighted graphs is faster
for large datasets (\> 100k cells).

``` {r}
#| fig-height: 6
#| fig-width: 6

pheatmap(alldata@graphs$CCA_nn[1:200,1:200],
         col=c("white","black"),border_color = "grey90", main = "KNN graph",
         legend = F,cluster_rows = F,cluster_cols = F,fontsize = 2)

pheatmap(alldata@graphs$CCA_snn[1:200,1:200],
         col=colorRampPalette(c("white","yellow","red"))(100), 
         border_color = "grey90", main = "SNN graph",
         legend = F,cluster_rows = F,cluster_cols = F,fontsize = 2)
```

### Clustering on a graph

Once the graph is built, we can now perform graph clustering. The
clustering is done respective to a resolution which can be interpreted
as how coarse you want your cluster to be. Higher resolution means
higher number of clusters.

In **Seurat**, the function `FindClusters()` will do a graph-based
clustering using "Louvain" algorithim by default (`algorithm = 1`). To
use the leiden algorithm, you need to set it to `algorithm = 4`. See
`?FindClusters` for additional options.

``` {r}
#| fig-height: 4
#| fig-width: 12
#| results: hide

# Clustering with louvain (algorithm 1)
for (res in c( 0.1 , 0.25 , .5 , 1 , 1.5 , 2 )){
  alldata <- FindClusters(alldata, graph.name = "CCA_snn", resolution = res , algorithm = 1)
}

# each time you run clustering, the data is stored in meta data columns:
# seurat_clusters - lastest results only
# CCA_snn_res.XX - for each different resolution you test.
```

``` {r}
#| fig-height: 4
#| fig-width: 12

plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.0.5")+ggtitle("louvain_0.5"),
  DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.1")+ggtitle("louvain_1"),
  DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.2")+ggtitle("louvain_2")
)
```

We can now use the `clustree` package to visualize how cells are
distributed between clusters depending on resolution.

``` {r}
#| fig-height: 8
#| fig-width: 8

suppressPackageStartupMessages(library(clustree))
clustree(alldata@meta.data, prefix = "CCA_snn_res.")
```

## K-means clustering

K-means is a generic clustering algorithm that has been used in many
application areas. In R, it can be applied via the kmeans function.
Typically, it is applied to a reduced dimension representation of the
expression data (most often PCA, because of the interpretability of the
low-dimensional distances). We need to define the number of clusters in
advance. Since the results depend on the initialization of the cluster
centers, it is typically recommended to run K-means with multiple
starting configurations (via the nstart argument).

``` {r}
#| fig-height: 4
#| fig-width: 12

for (k in c( 5 , 7 , 10 , 12 , 15 , 17 , 20)){
  alldata@meta.data[,paste0("kmeans_",k)] <- kmeans(x = alldata@reductions[["pca"]]@cell.embeddings, centers = k,nstart = 100)$cluster
}

plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "umap", group.by = "kmeans_5")+ggtitle("kmeans_5"),
  DimPlot(alldata, reduction = "umap", group.by = "kmeans_10")+ggtitle("kmeans_10"),
  DimPlot(alldata, reduction = "umap", group.by = "kmeans_15")+ggtitle("kmeans_15")
)
```

``` {r}
#| fig-height: 8
#| fig-width: 8

clustree(alldata@meta.data, prefix = "kmeans_")
```

## Hierarchical clustering

### Defining distance between cells

The base R `stats` package already contains a function `dist` that
calculates distances between all pairs of samples. Since we want to
compute distances between samples, rather than among genes, we need to
transpose the data before applying it to the `dist` function. This can
be done by simply adding the transpose function `t()` to the data. The
distance methods available in `dist` are: 'euclidean', 'maximum',
'manhattan', 'canberra', 'binary' or 'minkowski'.

``` {r}
d <- dist( alldata@reductions[["pca"]]@cell.embeddings, method="euclidean")
```

As you might have realized, correlation is not a method implemented in
the `dist()` function. However, we can create our own distances and
transform them to a distance object. We can first compute sample
correlations using the `cor` function.\
As you already know, correlation range from -1 to 1, where 1 indicates
that two samples are closest, -1 indicates that two samples are the
furthest and 0 is somewhat in between. This, however, creates a problem
in defining distances because a distance of 0 indicates that two samples
are closest, 1 indicates that two samples are the furthest and distance
of -1 is not meaningful. We thus need to transform the correlations to a
positive scale (a.k.a. **adjacency**):\
\[adj = `\frac{1- cor}{2}`{=tex}\]\
Once we transformed the correlations to a 0-1 scale, we can simply
convert it to a distance object using `as.dist` function. The
transformation does not need to have a maximum of 1, but it is more
intuitive to have it at 1, rather than at any other number.

``` {r}
#Compute sample correlations
sample_cor <- cor( Matrix::t(alldata@reductions[["pca"]]@cell.embeddings) )

#Transform the scale from correlations
sample_cor <- (1 - sample_cor) / 2

#Convert it to a distance object
d2 <- as.dist(sample_cor)
```

### Clustering cells

After having calculated the distances between samples calculated, we can
now proceed with the hierarchical clustering per-se. We will use the
function `hclust` for this purpose, in which we can simply run it with
the distance objects created above. The methods available are: 'ward.D',
'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or
'centroid'. It is possible to plot the dendrogram for all cells, but
this is very time consuming and we will omit for this tutorial.

``` {r}
#euclidean
h_euclidean <- hclust(d, method="ward.D2")

#correlation
h_correlation <- hclust(d2, method="ward.D2")
```

Once your dendrogram is created, the next step is to define which
samples belong to a particular cluster. After identifying the
dendrogram, we can now literally cut the tree at a fixed threshold (with
`cutree`) at different levels to define the clusters. We can either
define the number of clusters or decide on a height. We can simply try
different clustering levels.

``` {r}
#| fig-height: 8
#| fig-width: 13

#euclidean distance
alldata$hc_euclidean_5 <- cutree(h_euclidean,k = 5)
alldata$hc_euclidean_10 <- cutree(h_euclidean,k = 10)
alldata$hc_euclidean_15 <- cutree(h_euclidean,k = 15)

#correlation distance
alldata$hc_corelation_5 <- cutree(h_correlation,k = 5)
alldata$hc_corelation_10 <- cutree(h_correlation,k = 10)
alldata$hc_corelation_15 <- cutree(h_correlation,k = 15)


plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "umap", group.by = "hc_euclidean_5")+ggtitle("hc_euc_5"),
  DimPlot(alldata, reduction = "umap", group.by = "hc_euclidean_10")+ggtitle("hc_euc_10"),
  DimPlot(alldata, reduction = "umap", group.by = "hc_euclidean_15")+ggtitle("hc_euc_15"),

  DimPlot(alldata, reduction = "umap", group.by = "hc_corelation_5")+ggtitle("hc_cor_5"),
  DimPlot(alldata, reduction = "umap", group.by = "hc_corelation_10")+ggtitle("hc_cor_10"),
  DimPlot(alldata, reduction = "umap", group.by = "hc_corelation_15")+ggtitle("hc_cor_15")
)
```

Finally, lets save the clustered data for further analysis.

``` {r}
saveRDS(alldata,"data/results/covid_qc_dr_int_cl.rds")
```

<div>

> **Discuss**
>
> By now you should know how to plot different features onto your data.
> Take the QC metrics that were calculated in the first exercise, that
> should be stored in your data object, and plot it as violin plots per
> cluster using the clustering method of your choice. For example, plot
> number of UMIS, detected genes, percent mitochondrial reads. Then,
> check carefully if there is any bias in how your data is separated due
> to quality metrics. Could it be explained biologically, or could you
> have technical bias there?

</div>

## Session info

``` {r}
sessionInfo()
```
