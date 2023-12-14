---
description: Quality control of single cell RNA-Seq data. Inspection of
subtitle:  SCANPY TOOLKIT
title:  Quality Control
---

<div>

> **Note**
>
> Code chunks run Python commands unless it starts with `%%bash`, in
> which case, those chunks run shell commands.

</div>

## Get data

In this tutorial, we will run all tutorials with a set of 6 PBMC 10x
datasets from 3 covid-19 patients and 3 healthy controls, the samples
have been subsampled to 1500 cells per sample. They are part of the
github repo and if you have cloned the repo they should be available in
folder: `labs/data/covid_data_GSE149689`. Instructions on how to
download them can also be found in the Precourse material.

``` {python}
%%bash

if [ ! -d "data/raw" ] 
then
  # create a data directory.
  mkdir -p data/raw

  # first check if the files are there
  count=$(ls -l data/raw/*.h5 | grep -v ^d | wc -l )
  echo $count

  # if not 4 files, fetch the files from github.
  if (("$count" <  6)); then
    cd data/raw
    curl -O https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/new_dataset/labs/data/covid_data_GSE149689/sub/Normal_PBMC_13.h5
    curl -O https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/new_dataset/labs/data/covid_data_GSE149689/sub/Normal_PBMC_14.h5
    curl -O https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/new_dataset/labs/data/covid_data_GSE149689/sub/Normal_PBMC_5.h5
    curl -O https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/new_dataset/labs/data/covid_data_GSE149689/sub/nCoV_PBMC_15.h5
    curl -O https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/new_dataset/labs/data/covid_data_GSE149689/sub/nCoV_PBMC_17.h5
    curl -O https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/new_dataset/labs/data/covid_data_GSE149689/sub/nCoV_PBMC_1.h5
    cd ../..
  fi
fi

ls -lGa data/raw
```

With data in place, now we can start loading libraries we will use in
this tutorial.

``` {python}
import numpy as np
import pandas as pd
import scanpy as sc

# verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80)
```

We can first load the data individually by reading directly from HDF5
file format (.h5).

``` {python}
data_cov1 = sc.read_10x_h5('./data/raw/nCoV_PBMC_1.h5')
data_cov1.var_names_make_unique()
data_cov15 = sc.read_10x_h5('./data/raw/nCoV_PBMC_15.h5')
data_cov15.var_names_make_unique()
data_cov17 = sc.read_10x_h5('./data/raw/nCoV_PBMC_17.h5')
data_cov17.var_names_make_unique()
data_ctrl5 = sc.read_10x_h5('./data/raw/Normal_PBMC_5.h5')
data_ctrl5.var_names_make_unique()
data_ctrl13 = sc.read_10x_h5('./data/raw/Normal_PBMC_13.h5')
data_ctrl13.var_names_make_unique()
data_ctrl14 = sc.read_10x_h5('./data/raw/Normal_PBMC_14.h5')
data_ctrl14.var_names_make_unique()
```

## Collate

``` {python}
# add some metadata
data_cov1.obs['type']="Covid"
data_cov1.obs['sample']="covid_1"
data_cov15.obs['type']="Covid"
data_cov15.obs['sample']="covid_15"
data_cov17.obs['type']="Covid"
data_cov17.obs['sample']="covid_17"
data_ctrl5.obs['type']="Ctrl"
data_ctrl5.obs['sample']="ctrl_5"
data_ctrl13.obs['type']="Ctrl"
data_ctrl13.obs['sample']="ctrl_13"
data_ctrl14.obs['type']="Ctrl"
data_ctrl14.obs['sample']="ctrl_14"

# merge into one object.
adata = data_cov1.concatenate(data_cov15, data_cov17, data_ctrl5, data_ctrl13, data_ctrl14)

# and delete individual datasets to save space
del(data_cov1, data_cov15, data_cov17)
del(data_ctrl5, data_ctrl13, data_ctrl14)
```

You can print a summary of the datasets in the Scanpy object, or a
summary of the whole object.

``` {python}
print(adata.obs['sample'].value_counts())
adata
```

## Calculate QC

Having the data in a suitable format, we can start calculating some
quality metrics. We can for example calculate the percentage of
mitochondrial and ribosomal genes per cell and add to the metadata. The
proportion hemoglobin genes can give an indication of red blood cell
contamination. This will be helpful to visualize them across different
metadata parameteres (i.e. datasetID and chemistry version). There are
several ways of doing this. The QC metrics are finally added to the
metadata table.

Citing from Simple Single Cell workflows (Lun, McCarthy & Marioni,
2017): High proportions are indicative of poor-quality cells (Islam et
al. 2014; Ilicic et al. 2016), possibly because of loss of cytoplasmic
RNA from perforated cells. The reasoning is that mitochondria are larger
than individual transcript molecules and less likely to escape through
tears in the cell membrane.

First, let Scanpy calculate some general qc-stats for genes and cells
with the function `sc.pp.calculate_qc_metrics`, similar to
`calculateQCmetrics()` in Scater. It can also calculate proportion of
counts for specific gene populations, so first we need to define which
genes are mitochondrial, ribosomal and hemoglobin.

``` {python}
# mitochondrial genes
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
# ribosomal genes
adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
# hemoglobin genes.
adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)]"))

adata.var
```

``` {python}
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)
```

Now you can see that we have additional data in the metadata slot.

``` {python}
mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mt2'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

adata
```

## Plot QC

Now we can plot some of the QC variables as violin plots.

``` {python}
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo', 'pct_counts_hb'], jitter=0.4, groupby = 'sample', rotation= 45)
```

As you can see, there is quite some difference in quality for the 4
datasets, with for instance the covid_15 sample having fewer cells with
many detected genes and more mitochondrial content. As the ribosomal
proteins are highly expressed they will make up a larger proportion of
the transcriptional landscape when fewer of the lowly expressed genes
are detected. And we can plot the different QC-measures as scatter
plots.

``` {python}
#| fig-height: 5
#| fig-width: 5
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color="sample")
```

<div>

> **Discuss**
>
> Plot additional QC stats that we have calculated as scatter plots. How
> are the different measures correlated? Can you explain why?

</div>

## Filtering

### Detection-based filtering

A standard approach is to filter cells with low amount of reads as well
as genes that are present in at least a certain amount of cells. Here we
will only consider cells with at least 200 detected genes and genes need
to be expressed in at least 3 cells. Please note that those values are
highly dependent on the library preparation method used.

``` {python}
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(adata.n_obs, adata.n_vars)
```

Extremely high number of detected genes could indicate doublets.
However, depending on the cell type composition in your sample, you may
have cells with higher number of genes (and also higher counts) from one
cell type. In this case, we will run doublet prediction further down, so
we will skip this step now, but the code below is an example of how it
can be run:

``` {python}
# skip for now as we are doing doublet prediction
#keep_v2 = (adata.obs['n_genes_by_counts'] < 2000) & (adata.obs['n_genes_by_counts'] > 500) & (adata.obs['lib_prep'] == 'v2')
#print(sum(keep_v2))

# filter for gene detection for v3
#keep_v3 = (adata.obs['n_genes_by_counts'] < 4100) & (adata.obs['n_genes_by_counts'] > 1000) & (adata.obs['lib_prep'] != 'v2')
#print(sum(keep_v3))

# keep both sets of cells
#keep = (keep_v2) | (keep_v3)
#print(sum(keep))
#adata = adata[keep, :]

#print("Remaining cells %d"%adata.n_obs)
```

Additionally, we can also see which genes contribute the most to such
reads. We can for instance plot the percentage of counts per gene.

``` {python}
#| fig-height: 6
#| fig-width: 6
sc.pl.highest_expr_genes(adata, n_top=20)
```

As you can see, MALAT1 constitutes up to 30% of the UMIs from a single
cell and the other top genes are mitochondrial and ribosomal genes. It
is quite common that nuclear lincRNAs have correlation with quality and
mitochondrial reads, so high detection of MALAT1 may be a technical
issue. Let us assemble some information about such genes, which are
important for quality control and downstream filtering.

### Mito/Ribo filtering

We also have quite a lot of cells with high proportion of mitochondrial
and low proportion of ribosomal reads. It could be wise to remove those
cells, if we have enough cells left after filtering. Another option
would be to either remove all mitochondrial reads from the dataset and
hope that the remaining genes still have enough biological signal. A
third option would be to just regress out the `percent_mito` variable
during scaling. In this case we had as much as 99.7% mitochondrial reads
in some of the cells, so it is quite unlikely that there is much cell
type signature left in those. Looking at the plots, make reasonable
decisions on where to draw the cutoff. In this case, the bulk of the
cells are below 20% mitochondrial reads and that will be used as a
cutoff. We will also remove cells with less than 5% ribosomal reads.

``` {python}
# filter for percent mito
adata = adata[adata.obs['pct_counts_mt'] < 20, :]

# filter for percent ribo > 0.05
adata = adata[adata.obs['pct_counts_ribo'] > 5, :]

print("Remaining cells %d"%adata.n_obs)
```

As you can see, a large proportion of sample covid_15 is filtered out.
Also, there is still quite a lot of variation in `percent_mito`, so it
will have to be dealt with in the data analysis step. We can also notice
that the `percent_ribo` are also highly variable, but that is expected
since different cell types have different proportions of ribosomal
content, according to their function.

### Plot filtered QC

Lets plot the same QC-stats another time.

``` {python}
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo', 'pct_counts_hb'], jitter=0.4, groupby = 'sample', rotation = 45)
```

### Filter genes

As the level of expression of mitochondrial and MALAT1 genes are judged
as mainly technical, it can be wise to remove them from the dataset
before any further analysis.

``` {python}
malat1 = adata.var_names.str.startswith('MALAT1')
# we need to redefine the mito_genes since they were first 
# calculated on the full object before removing low expressed genes.
mito_genes = adata.var_names.str.startswith('MT-')
hb_genes = adata.var_names.str.contains('^HB[^(P)]')

remove = np.add(mito_genes, malat1)
remove = np.add(remove, hb_genes)
keep = np.invert(remove)

adata = adata[:,keep]

print(adata.n_obs, adata.n_vars)
```

## Sample sex

When working with human or animal samples, you should ideally constrain
you experiments to a single sex to avoid including sex bias in the
conclusions. However this may not always be possible. By looking at
reads from chromosomeY (males) and XIST (X-inactive specific transcript)
expression (mainly female) it is quite easy to determine per sample
which sex it is. It can also bee a good way to detect if there has been
any sample mixups, if the sample metadata sex does not agree with the
computational predictions.

To get choromosome information for all genes, you should ideally parse
the information from the gtf file that you used in the mapping pipeline
as it has the exact same annotation version/gene naming. However, it may
not always be available, as in this case where we have downloaded public
data. Hence, we will use biomart to fetch chromosome information.

``` {python}
# requires pybiomart
annot = sc.queries.biomart_annotations("hsapiens", ["ensembl_gene_id", "external_gene_name", "start_position", "end_position", "chromosome_name"], ).set_index("external_gene_name")
# adata.var[annot.columns] = annot
```

Now that we have the chromosome information, we can calculate per cell
the proportion of reads that comes from chromosome Y.

``` {python}
chrY_genes = adata.var_names.intersection(annot.index[annot.chromosome_name == "Y"])
chrY_genes

adata.obs['percent_chrY'] = np.sum(
    adata[:, chrY_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1 * 100
```

Then plot XIST expression vs chrY proportion. As you can see, the
samples are clearly on either side, even if some cells do not have
detection of either.

``` {python}
#| fig-height: 5
#| fig-width: 5

# color inputs must be from either .obs or .var, so add in XIST expression to obs.
adata.obs["XIST-counts"] = adata.X[:,adata.var_names.str.match('XIST')].toarray()

sc.pl.scatter(adata, x='XIST-counts', y='percent_chrY', color="sample")
```

Plot as violins.

``` {python}
#| fig-height: 5
#| fig-width: 10

sc.pl.violin(adata, ["XIST-counts", "percent_chrY"], jitter=0.4, groupby = 'sample', rotation= 45)
```

Here, we can see clearly that we have two males and 4 females, can you
see which samples they are? Do you think this will cause any problems
for downstream analysis? Discuss with your group: what would be the best
way to deal with this type of sex bias?

## Cell cycle state

We here perform cell cycle scoring. To score a gene list, the algorithm
calculates the difference of mean expression of the given list and the
mean expression of reference genes. To build the reference, the function
randomly chooses a bunch of genes matching the distribution of the
expression of the given list. Cell cycle scoring adds three slots in
data, a score for S phase, a score for G2M phase and the predicted cell
cycle phase.

First read the file with cell cycle genes, from Regev lab and split into
S and G2M phase genes. Cell cycle genes were retrieved from the
scanpy_usage github site via web browser at [RegevLab Github
repo](https://github.com/theislab/scanpy_usage/blob/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt).

``` {python}
%%bash
if [ ! -f data/regev_lab_cell_cycle_genes.txt ]; then curl -o data/regev_lab_cell_cycle_genes.txt https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt; fi
```

``` {python}
cell_cycle_genes = [x.strip() for x in open('./data/regev_lab_cell_cycle_genes.txt')]
print(len(cell_cycle_genes))

# Split into 2 lists
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
print(len(cell_cycle_genes))
```

Before running cell cycle we have to normalize the data. In the scanpy
object, the data slot will be overwritten with the normalized data. So
first, save the raw data into the slot `raw`. Then run normalization,
log transformation and scale the data.

``` {python}
# save normalized counts in raw slot.
adata.raw = adata

# normalize to depth 10 000
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

# logaritmize
sc.pp.log1p(adata)

# scale
sc.pp.scale(adata)
```

We here perform cell cycle scoring. The function is actually a wrapper
to sc.tl.score_gene_list, which is launched twice, to score separately S
and G2M phases. Both sc.tl.score_gene_list and
sc.tl.score_cell_cycle_genes are a port from Seurat and are supposed to
work in a very similar way. To score a gene list, the algorithm
calculates the difference of mean expression of the given list and the
mean expression of reference genes. To build the reference, the function
randomly chooses a bunch of genes matching the distribution of the
expression of the given list. Cell cycle scoring adds three slots in
data, a score for S phase, a score for G2M phase and the predicted cell
cycle phase.

``` {python}
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
```

We can now plot a violin plot for the cell cycle scores as well.

``` {python}
#| fig-height: 5
#| fig-width: 10

sc.pl.violin(adata, ['S_score', 'G2M_score'], jitter=0.4, groupby = 'sample', rotation=45)
```

In this case it looks like we only have a few cycling cells in the
datasets.

## Predict doublets

Doublets/Multiples of cells in the same well/droplet is a common issue
in scRNAseq protocols. Especially in droplet-based methods with
overloading of cells. In a typical 10x experiment the proportion of
doublets is linearly dependent on the amount of loaded cells. As
indicated from the Chromium user guide, doublet rates are about as
follows:\
![](../figs/10x_doublet_rate.png)\
Most doublet detectors simulates doublets by merging cell counts and
predicts doublets as cells that have similar embeddings as the simulated
doublets. Most such packages need an assumption about the
number/proportion of expected doublets in the dataset. The data you are
using is subsampled, but the original datasets contained about 5 000
cells per sample, hence we can assume that they loaded about 9 000 cells
and should have a doublet rate at about 4%.

<div>

> **Caution**
>
> Ideally doublet prediction should be run on each sample separately,
> especially if your different samples have different proportions of
> cell types. In this case, the data is subsampled so we have very few
> cells per sample and all samples are sorted PBMCs so it is okay to run
> them together.

</div>

For doublet detection, we will use the package `Scrublet`, so first we
need to get the raw counts from `adata.raw.X` and run scrublet with that
matrix. Then we add in the doublet prediction info into our anndata
object.

Doublet prediction should be run for each dataset separately, so first
we need to split the adata object into 6 separate objects, one per
sample and then run scrublet on each of them.

``` {python}
import scrublet as scr

# split per batch into new objects.
batches = adata.obs['sample'].cat.categories.tolist()
alldata = {}
for batch in batches:
    tmp = adata[adata.obs['sample'] == batch,]
    print(batch, ":", tmp.shape[0], " cells")
    scrub = scr.Scrublet(tmp.raw.X)
    out = scrub.scrub_doublets(verbose=False, n_prin_comps = 20)
    alldata[batch] = pd.DataFrame({'doublet_score':out[0],'predicted_doublets':out[1]},index = tmp.obs.index)
    print(alldata[batch].predicted_doublets.sum(), " predicted_doublets")
```

``` {python}
# add predictions to the adata object.
scrub_pred = pd.concat(alldata.values())
adata.obs['doublet_scores'] = scrub_pred['doublet_score'] 
adata.obs['predicted_doublets'] = scrub_pred['predicted_doublets'] 

sum(adata.obs['predicted_doublets'])
```

We should expect that two cells have more detected genes than a single
cell, lets check if our predicted doublets also have more detected genes
in general.

``` {python}
#| fig-height: 5
#| fig-width: 5

# add in column with singlet/doublet instead of True/Fals
%matplotlib inline

adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)
sc.pl.violin(adata, 'n_genes_by_counts', jitter=0.4, groupby = 'doublet_info', rotation=45)
```

Now, lets run PCA and UMAP and plot doublet scores onto UMAP to check
the doublet predictions.

``` {python}
#| fig-height: 4
#| fig-width: 12

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['doublet_scores','doublet_info','sample'])
```

Now, lets remove all predicted doublets from our data.

``` {python}
# also revert back to the raw counts as the main matrix in adata
adata = adata.raw.to_adata() 

adata = adata[adata.obs['doublet_info'] == 'False',:]
print(adata.shape)
```

## Save data

Finally, lets save the QC-filtered data for further analysis. Create
output directory `results` and save data to that folder. This will be
used in downstream labs.

``` {python}
import os
if not os.path.exists('data/results/'):
    os.makedirs('data/results/', exist_ok=True)
if not os.path.exists('data/results/scanpy_qc_filtered_covid.h5ad'):
    adata.write_h5ad('data/results/scanpy_qc_filtered_covid.h5ad')
```

## Session info

``` {python}
sc.logging.print_versions()
```
