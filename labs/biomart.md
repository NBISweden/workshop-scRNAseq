biomaRt
=======

Author: Åsa Björklund

Here is some example code on how to translate between gene symbols and ensembl gene ids using the biomaRt package. For more details on the package, have a look at: <https://bioconductor.org/packages/release/bioc/html/biomaRt.html>

All data you need is available in the course uppmax folder with subfolder: `scrnaseq_course/data/ILC/`

#### Select dataset

Load the biomaRt package and select which mart and dataset to use.

``` r
suppressMessages(require(biomaRt))

# select which mart to use, in this case ensembl
mart <- useMart("ensembl")

# To see what datasets exits you can run: listDatasets
head(listDatasets(mart))
```

    ##                   dataset                              description
    ## 1 amexicanus_gene_ensembl Cave fish genes (Astyanax_mexicanus-2.0)
    ## 2  mmusculus_gene_ensembl                  Mouse genes (GRCm38.p6)
    ## 3 nbrichardi_gene_ensembl       Lyretail cichlid genes (NeoBri1.0)
    ## 4  pmexicana_gene_ensembl    Shortfin molly genes (P_mexicana-1.0)
    ## 5 lchalumnae_gene_ensembl               Coelacanth genes (LatCha1)
    ## 6      mfuro_gene_ensembl              Ferret genes (MusPutFur1.0)
    ##                  version
    ## 1 Astyanax_mexicanus-2.0
    ## 2              GRCm38.p6
    ## 3              NeoBri1.0
    ## 4         P_mexicana-1.0
    ## 5                LatCha1
    ## 6           MusPutFur1.0

``` r
# in this case we use hsapiens_gene_ensembl
mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
```

#### Search based on Ensembl ID

Here we will fetch gene\_id, gene\_name, description, biotype and chromosome position for all ensembl\_ids that we have in the expression matrix.

``` r
# to find out what attributes there are in the Dataset, use listAttributes
head(listAttributes(mart))
```

    ##                            name                  description         page
    ## 1               ensembl_gene_id               Gene stable ID feature_page
    ## 2       ensembl_gene_id_version       Gene stable ID version feature_page
    ## 3         ensembl_transcript_id         Transcript stable ID feature_page
    ## 4 ensembl_transcript_id_version Transcript stable ID version feature_page
    ## 5            ensembl_peptide_id            Protein stable ID feature_page
    ## 6    ensembl_peptide_id_version    Protein stable ID version feature_page

``` r
# read in expression matrix to get the genes we want to translate
R <- read.table("data/ILC/ensembl_rpkmvalues_ILC.csv",sep=",",header=T)

# getBM function fetches attributes from the database with specified names. 
# with filters parameter you define which attribute you want to filter on
# with values, you define which entries you want to fetch, leave empty to fetch all entries.
# with attributes, you define what attributes you want to fetch

genes.table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", 
              "external_gene_name", "description","gene_biotype", "chromosome_name","start_position"), 
              values= rownames(R), mart= mart) 
```

    ## Warning in strptime(x, fmt, tz = "GMT"): unknown timezone 'zone/tz/2018i.
    ## 1.0/zoneinfo/Europe/Stockholm'

``` r
head(genes.table)
```

    ##   ensembl_gene_id external_gene_name
    ## 1 ENSG00000000003             TSPAN6
    ## 2 ENSG00000000005               TNMD
    ## 3 ENSG00000000419               DPM1
    ## 4 ENSG00000000457              SCYL3
    ## 5 ENSG00000000460           C1orf112
    ## 6 ENSG00000000938                FGR
    ##                                                                                      description
    ## 1                                              tetraspanin 6 [Source:HGNC Symbol;Acc:HGNC:11858]
    ## 2                                                tenomodulin [Source:HGNC Symbol;Acc:HGNC:17757]
    ## 3 dolichyl-phosphate mannosyltransferase subunit 1, catalytic [Source:HGNC Symbol;Acc:HGNC:3005]
    ## 4                                   SCY1 like pseudokinase 3 [Source:HGNC Symbol;Acc:HGNC:19285]
    ## 5                        chromosome 1 open reading frame 112 [Source:HGNC Symbol;Acc:HGNC:25565]
    ## 6              FGR proto-oncogene, Src family tyrosine kinase [Source:HGNC Symbol;Acc:HGNC:3697]
    ##     gene_biotype chromosome_name start_position
    ## 1 protein_coding               X      100627109
    ## 2 protein_coding               X      100584802
    ## 3 protein_coding              20       50934867
    ## 4 protein_coding               1      169849631
    ## 5 protein_coding               1      169662007
    ## 6 protein_coding               1       27612064

``` r
# write to a file for later use
write.table(genes.table, file="data/ILC/gene_name_translation_biotype.tab",sep="\t")
```

#### Fetch Ensembl ID based on gene names

You can do the opposite if you have gene names and want Ensembl IDs.

``` r
# now we want to get all ensembl IDs for the genes in genes.table$external_gene_name
genes.table2 <- getBM(filters= "external_gene_name", attributes= c("ensembl_gene_id", 
          "external_gene_name", "description","gene_biotype"), 
          values= genes.table$external_gene_name, mart= mart)

# Keep in mind, you may get multiple ensembl IDs translated to the same gene name, 
# so the number of entries will be different.
dim(genes.table)
```

    ## [1] 48909     6

``` r
dim(genes.table2)
```

    ## [1] 69822     4

Also, keep in mind that if you are working with an older version of Ensembl, some Ensembl IDs may be obsolete and not have any translation, so those will require some manual searching to annotate with gene names.

#### Fetch Gene Ontology annotations

You can also use biomaRt to fetch gene ontology annotations and a bunch of other attributes that you can find in the database. Here is an example for fetching GO-terms, that may be useful for running Pagoda if you are using your own dataset.

``` r
go.table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name", "go_id","name_1006", "namespace_1003"), values= rownames(R), mart= mart)

head(go.table)
```

    ##   ensembl_gene_id external_gene_name      go_id
    ## 1 ENSG00000000003             TSPAN6 GO:0016020
    ## 2 ENSG00000000003             TSPAN6 GO:0016021
    ## 3 ENSG00000000003             TSPAN6 GO:0007166
    ## 4 ENSG00000000003             TSPAN6 GO:0005887
    ## 5 ENSG00000000003             TSPAN6 GO:0005515
    ## 6 ENSG00000000003             TSPAN6 GO:0070062
    ##                                 name_1006     namespace_1003
    ## 1                                membrane cellular_component
    ## 2          integral component of membrane cellular_component
    ## 3 cell surface receptor signaling pathway biological_process
    ## 4   integral component of plasma membrane cellular_component
    ## 5                         protein binding molecular_function
    ## 6                   extracellular exosome cellular_component

``` r
# If you want to create a list with all genes as keys, and a vector of go-terms as values
gene2go <- split(go.table$go_id, go.table$ensembl_gene_id)
head(gene2go)
```

    ## $ENSG00000000003
    ##  [1] "GO:0016020" "GO:0016021" "GO:0007166" "GO:0005887" "GO:0005515"
    ##  [6] "GO:0070062" "GO:0043123" "GO:0039532" "GO:1901223" ""          
    ## 
    ## $ENSG00000000005
    ##  [1] "GO:0005634" "GO:0005737" "GO:0016020" "GO:0016021" "GO:0005515"
    ##  [6] "GO:0016525" "GO:0005635" "GO:0001886" "GO:0001937" "GO:0030948"
    ## [11] "GO:0035990" "GO:0071773" ""          
    ## 
    ## $ENSG00000000419
    ##  [1] "GO:0005634" "GO:0016020" "GO:0016740" "GO:0005783" "GO:0016757"
    ##  [6] "GO:0005515" "GO:0006486" "GO:0006506" "GO:0005789" "GO:0035269"
    ## [11] "GO:0004169" "GO:0004582" "GO:0019348" "GO:0035268" "GO:0033185"
    ## [16] ""           "GO:0097502" "GO:0005537" "GO:0043178" "GO:0019673"
    ## [21] "GO:0043231"
    ## 
    ## $ENSG00000000457
    ##  [1] "GO:0005737" "GO:0005524" "GO:0006468" "GO:0005794" "GO:0042995"
    ##  [6] "GO:0005515" "GO:0030027" "GO:0016477" "GO:0004672" ""          
    ## 
    ## $ENSG00000000460
    ## [1] ""
    ## 
    ## $ENSG00000000938
    ##  [1] "GO:0005737" "GO:0016020" "GO:0000166" "GO:0004672" "GO:0005524"
    ##  [6] "GO:0006468" "GO:0016740" "GO:0016301" "GO:0016310" "GO:0005743"
    ## [11] "GO:0004713" "GO:0005886" "GO:0005829" "GO:0005739" "GO:0005856"
    ## [16] "GO:0042995" "GO:0005515" "GO:0016235" "GO:0004715" "GO:0007229"
    ## [21] "GO:0002376" "GO:0045087" "GO:0038096" "GO:0005576" "GO:0030154"
    ## [26] "GO:0007169" "GO:0005758" "GO:0031234" "GO:0070062" "GO:0019901"
    ## [31] "GO:0032587" "GO:0030335" "GO:0050830" "GO:0034987" "GO:0043312"
    ## [36] "GO:0008360" "GO:0034774" "GO:0046777" "GO:0001784" "GO:0005102"
    ## [41] "GO:0009615" "GO:0045859" "GO:0015629" "GO:0014068" "GO:0050764"
    ## [46] "GO:0018108" "GO:0038083" "GO:0043552" "GO:0050715" "GO:0002768"
    ## [51] "GO:0034988" "GO:0043306" "GO:0045088" ""

``` r
# To do the opposite, go-terms as keys with a vector of genes with that go-term
go2gene <- split(go.table$ensembl_gene_id, go.table$go_id)
```

#### Select only Biological Process

If you want to select only Biological Process, the entry "namespace\_1003" defines the type of GO-term, so you can filter on that as well.

``` r
go.tableBP <- go.table[go.table$namespace_1003=="biological_process",]

# if you want more informative names for the go-terms, merge GO-id with name
go.name <- paste(go.tableBP$go_id,go.tableBP$name_1006,sep=";")
go.tableBP$go.name <- go.name

#make a list with GO-name to gene IDs
goBP2gene <- split(go.tableBP$ensembl_gene_id, go.tableBP$go.name)

# save to file
save(goBP2gene, file="data/ILC/GO_BP_annotations.Rdata")
```

##### Session info

``` r
sessionInfo()
```

    ## R version 3.4.1 (2017-06-30)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS  10.14
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] biomaRt_2.34.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.15         AnnotationDbi_1.40.0 knitr_1.19          
    ##  [4] magrittr_1.5         progress_1.1.2       IRanges_2.12.0      
    ##  [7] BiocGenerics_0.24.0  bit_1.1-12           R6_2.2.2            
    ## [10] rlang_0.3.0.1        httr_1.3.1           stringr_1.3.1       
    ## [13] blob_1.1.0           tools_3.4.1          parallel_3.4.1      
    ## [16] Biobase_2.38.0       DBI_0.7              htmltools_0.3.6     
    ## [19] assertthat_0.2.0     yaml_2.1.16          bit64_0.9-7         
    ## [22] rprojroot_1.3-2      digest_0.6.15        tibble_1.4.2        
    ## [25] S4Vectors_0.16.0     bitops_1.0-6         curl_3.1            
    ## [28] RCurl_1.95-4.10      memoise_1.1.0        evaluate_0.10.1     
    ## [31] RSQLite_2.0          rmarkdown_1.8        stringi_1.2.4       
    ## [34] pillar_1.1.0         compiler_3.4.1       prettyunits_1.0.2   
    ## [37] backports_1.1.2      stats4_3.4.1         XML_3.98-1.9
