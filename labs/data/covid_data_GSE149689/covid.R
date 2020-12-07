# download to folder full:
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149689/matrix/GSE149689_series_matrix.txt.gz
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149689/suppl/GSE149689_barcodes.tsv.gz
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149689/suppl/GSE149689_features.tsv.gz
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149689/suppl/GSE149689_matrix.mtx.gz


#GSE149689_series_matrix
library(Matrix)
library(hdf5r)
library(rhdf5)

sm <- as(Matrix::readMM("full/GSE149689_matrix.mtx.gz"),Class = "dgCMatrix")
sm@Dimnames[[1]] <- as.character(read.delim("full/GSE149689_features.tsv.gz",header = F)[,2])
sm@Dimnames[[2]] <- as.character(read.delim("full/GSE149689_barcodes.tsv.gz",header = F)[,1])
dim(sm)



#see sample metadata
meta <- read.delim("full/GSE149689_series_matrix.txt.gz",skip = 54)
t(meta[c(7,8,9,10,11),])

#check how many cells per sample.
t = table(sub(".*-","",colnames(sm)))[as.character(1:length(table(sub(".*-","",colnames(sm)))) )]
cbind(t,colnames(meta)[2:21],t(unname(meta[10,2:21])))

# 1  "6455" "Sample.1_nCoV.1.scRNA.seq..SW103."    "subject status: severe COVID-19 patient"              
# 2  "7731" "Sample.2_nCoV.2.scRNA.seq..SW104."    "subject status: mild COVID-19 patient"                
# 3  "5851" "Sample.3_Flu.1.scRNA.seq..SW105."     "age: 68"                                              
# 4  "1724" "Sample.4_Flu.2.scRNA.seq..SW106."     "age: 75"                                              
# 5  "6426" "Sample.5_Normal.1.scRNA.seq..SW107."  "age: 63"                                              
# 6  "3516" "Sample.6_Flu.3.scRNA.seq..SW108."     "age: 70"                                              
# 7  "1455" "Sample.7_Flu.4.scRNA.seq..SW109."     "age: 56"                                              
# 8  "2001" "Sample.8_Flu.5.scRNA.seq..SW110."     "age: 78"                                              
# 9  "538"  "Sample.9_nCoV.3.scRNA.seq..SW111."    "subject status: severe COVID-19 patient"              
# 10 "1167" "Sample.10_nCoV.4.scRNA.seq..SW112."   "subject status: severe COVID-19 patient"              
# 11 "6398" "Sample.11_nCoV.5.scRNA.seq..SW113."   "subject status: mild COVID-19 patient"                
# 12 "6283" "Sample.12_nCoV.6.scRNA.seq..SW114."   "subject status: mild COVID-19 patient"                
# 13 "6156" "Sample.13_Normal.2.scRNA.seq..SW115." "age: 54"                                              
# 14 "6574" "Sample.14_Normal.3.scRNA.seq..SW116." "age: 67"                                              
# 15 "2349" "Sample.15_nCoV.7.scRNA.seq..SW117."   "subject status: severe COVID-19 patient"              
# 16 "4371" "Sample.16_nCoV.8.scRNA.seq..SW118."   "subject status: severe COVID-19 patient"              
# 17 "1755" "Sample.17_nCoV.9.scRNA.seq..SW119."   "subject status: severe COVID-19 patient"              
# 18 "3984" "Sample.18_nCoV.10.scRNA.seq..SW120."  "subject status: mild COVID-19 patient"                
# 19 "5542" "Sample.19_Normal.4.scRNA.seq..SW121." "age: 63"                                              
# 20 "4868" "Sample.20_nCoV.11.scRNA.seq..SW122."  "subject status: Asymptomatic case of COVID-19 patient"


# select severe covid samples with many cells. 15,16 probably best.

# filter for at least 200 umis per sample and recount, have cells with only 15 umis...
nC = colSums(sm)
range(nC)
dim(sm)
sm = sm[,nC>200]
dim(sm)
# removes 708 cells. Mainly from samples 8,10,17.

t2 = table(sub(".*-","",colnames(sm)))[as.character(1:length(table(sub(".*-","",colnames(sm)))) )]
cbind(t2,colnames(meta)[2:21],t(unname(meta[10,2:21])))

# 1  "6455" "Sample.1_nCoV.1.scRNA.seq..SW103."    "subject status: severe COVID-19 patient"              
# 2  "7731" "Sample.2_nCoV.2.scRNA.seq..SW104."    "subject status: mild COVID-19 patient"                
# 3  "5835" "Sample.3_Flu.1.scRNA.seq..SW105."     "age: 68"                                              
# 4  "1656" "Sample.4_Flu.2.scRNA.seq..SW106."     "age: 75"                                              
# 5  "6426" "Sample.5_Normal.1.scRNA.seq..SW107."  "age: 63"                                              
# 6  "3501" "Sample.6_Flu.3.scRNA.seq..SW108."     "age: 70"                                              
# 7  "1382" "Sample.7_Flu.4.scRNA.seq..SW109."     "age: 56"                                              
# 8  "1804" "Sample.8_Flu.5.scRNA.seq..SW110."     "age: 78"                                              
# 9  "481"  "Sample.9_nCoV.3.scRNA.seq..SW111."    "subject status: severe COVID-19 patient"              
# 10 "1033" "Sample.10_nCoV.4.scRNA.seq..SW112."   "subject status: severe COVID-19 patient"              
# 11 "6398" "Sample.11_nCoV.5.scRNA.seq..SW113."   "subject status: mild COVID-19 patient"                
# 12 "6283" "Sample.12_nCoV.6.scRNA.seq..SW114."   "subject status: mild COVID-19 patient"                
# 13 "6156" "Sample.13_Normal.2.scRNA.seq..SW115." "age: 54"                                              
# 14 "6574" "Sample.14_Normal.3.scRNA.seq..SW116." "age: 67"                                              
# 15 "2290" "Sample.15_nCoV.7.scRNA.seq..SW117."   "subject status: severe COVID-19 patient"              
# 16 "4371" "Sample.16_nCoV.8.scRNA.seq..SW118."   "subject status: severe COVID-19 patient"              
# 17 "1666" "Sample.17_nCoV.9.scRNA.seq..SW119."   "subject status: severe COVID-19 patient"              
# 18 "3984" "Sample.18_nCoV.10.scRNA.seq..SW120."  "subject status: mild COVID-19 patient"                
# 19 "5542" "Sample.19_Normal.4.scRNA.seq..SW121." "age: 63"                                              
# 20 "4868" "Sample.20_nCoV.11.scRNA.seq..SW122."  "subject status: Asymptomatic case of COVID-19 patient"


# select all severe and normal samples and create matrices. Skip sample 9,10, less than 1500 cells. 
samples_use <- c(c(1,15,16,17),c(5,13,14,19))
sum(table(sub(".*-","",colnames(sm)))[as.character(samples_use)])



#subset dataset
# sm2 <- sm[,grepl(paste("-",samples_use,"$",sep = "",collapse = "|"),colnames(sm))]
# dim(sm2)
sel <- unlist(lapply(samples_use,function(x){
  set.seed(1);x <- sample(size = 1500,grep(paste0("-",x,"$"),colnames(sm),value = T) )
}))
sm2 <- sm[,sel]
table(sub(".*-","",colnames(sm2)))[as.character(samples_use)]

dim(sm2)


# need to add in gene id column as well.
feats = read.delim("full/GSE149689_features.tsv.gz",header = F)



for(i in c(paste0("nCoV_PBMC_",c(1,15,16,17)), paste0("Normal_PBMC_",c(5,13,14,19)) )) {
  message(paste0("PROCESSING SAMPLE:    ",i) )
  spn <- sub(".*_","",i)
  fn <- paste0("sub/",i,".h5")
  group <- grep(paste0("-",spn,"$"),colnames(sm2),value = T)
  sm3 <-  sm2[,group]
  dim(sm3)
  
  # file.remove(fn)
  rhdf5::h5createFile(fn)
  rhdf5::h5createGroup(fn,"matrix")
  
  rhdf5::h5write(sm3@Dimnames[[2]],fn,"matrix/barcodes")
  rhdf5::h5write(sm3@x,fn,"matrix/data")
  rhdf5::h5write(sm3@i,fn,"matrix/indices")
  rhdf5::h5write(sm3@p,fn,"matrix/indptr")
  rhdf5::h5write(sm3@Dim,fn,"matrix/shape")
  
  rhdf5::h5createGroup(fn,"matrix/features")
  rhdf5::h5write(sm3@Dimnames[[1]]
          ,fn,"matrix/features/name")
  rhdf5::h5write(sm3@Dimnames[[1]]
          ,fn,"matrix/features/_all_tag_keys")
  rhdf5::h5write(feats[,3],
          fn,"matrix/features/feature_type")
  rhdf5::h5write(feats[,1],
          fn,"matrix/features/id")
  rhdf5::h5write(rep("GRCh38",nrow(sm3))
          ,fn,"matrix/features/genome")
  
  rhdf5::h5ls(fn)
  
  nd <- Seurat::Read10X_h5(fn)
  print(sum(!nd == sm3))
  message("\n\n")
}


# check quality of all

library(Seurat)
sid = sub(".*-","",colnames(sm2))
m = data.frame(sample = paste0("sm",sid), type = ifelse(sid %in% c(5,13,14,19), "Normal","Covid"))
rownames(m) = colnames(sm2)
m$sample2 = paste(m$type, m$sample, sep="_")
sdata = CreateSeuratObject(sm2, meta.data = m)

pdf("sample_qualities.pdf")
VlnPlot(sdata, features = c("nFeature_RNA", "nCount_RNA"), group.by = "sample2", pt.size = 0)
dev.off()


# select covid 1,15,17 and normal 13,14,5




