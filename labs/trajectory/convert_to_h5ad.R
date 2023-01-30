# convert Seurat data to h5ad in R

library(Seurat)
library(SeuratDisk)

tdata = readRDS("./../data/bone_marrow/trajectory_seurat_filtered.rds")

#OBS! Scanpy will consider columns with numbers as factors as numerical, need to convert them to strings. For instance the cluster numbers.
tdata$clusters = as.character(tdata$clusters)
tdata$metadata_clusters = as.character(tdata$metadata_clusters)


tmpfile = "./../data/bone_marrow/trajectory_scanpy_filtered.h5Seurat"

SaveH5Seurat(tdata, filename = tmpfile)
Convert(tmpfile, dest = "h5ad", overwrite = T)

file.remove(tmpfile)

