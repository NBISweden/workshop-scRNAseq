# Convert labs to include tag text
# render as new ipynb and html for scanpy
# render as new Rmd and html for scater/seurat

# OBS! Requires jupyter_contrib_nbextensions to be installed.
# OBS! scanpy_07_spatial.ipynb requires a different conda env, has conflicts with nbconvert version.

# to run conversion of all labs just run
# Rscript convert_all_labs.R -f all -t all


# input is either --file all, of --file filename
# flag -t should be scater, scran, seurat or all

library(optparse)
library(RJSONIO)

#parse input
option_list <- list( 
  make_option(c("-f", "--file"), type = "character", default="all",
              help="Script to convert, if or leave as all to run for all scripts [all]"),
  make_option(c("-t", "--type"), type = "character", default="scanpy",
              help="Which type of lab, options (scanpy, scater, seurat), if or leave as all to run all of them [all]")  
)
opt <- parse_args(OptionParser(option_list=option_list))


#define path to this script 
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
setwd(script_path)

# for testing...
#opt <- list()
#opt$file <- "scanpy_04_clustering.ipynb"
#script_path <- "."

#opt <- list()
#opt$file <- "scater_01_qc.Rmd"
#opt$type <- "scater"
#script_path <- "."

# read the ref text
lab_text <- readLines(paste0(script_path,"/knit_labs.Rmd"))

######################################
# Scanpy
######################################

if (opt$type == "all" | opt$type == "scanpy"){
  
  # create output dir
  outdir <- "compiled/scanpy"
  dir.create(outdir,recursive = T, showWarnings = F)

  # check for input files
  indir <- "scanpy"
  if (opt$file == "all"){
    scripts <- dir(indir, pattern = "*.ipynb")
  }else {
    scripts <- opt$file
    if (!file.exists(file.path(indir,scripts))){
      stop(sprintf("File %s does not exist in dir %s",scripts,indir))
    }
  }



# function to replace text in cells
replace.tags <- function(in_text, ref_text){
  ipy_m <- regexpr("[#].*[_].*[:]",in_text,perl = T)
  ipy_tags <- regmatches(in_text, ipy_m)
  ref_tags <- sub(":.*",":",grep("[#].*[_].*[:]",ref_text,value = T))
  
  matched_tags <- ipy_tags %in% ref_tags
  
  if(any(!matched_tags)){
    for (t in ipy_tags[!matched_tags]){
      print(sprintf("No matched tag for %s",t))
    }
  }
  t_lab <- as.list(in_text)
  for( tag in ipy_tags[matched_tags] ){
    j <- grep(tag,t_lab)
    replace <- lab_text[grepl(tag,ref_text)]
    replace <- gsub(tag,"",replace)
    # add \n to each line
    replace <- paste0(replace, "\n")
    t_lab[[j]] <- replace
  }
  return(unlist(t_lab))
}

for (file in scripts){
  infile <- file.path(indir,file)
  print(sprintf("Parsing data from %s",infile))
  json <- RJSONIO::fromJSON(infile)

  #json_temp <- json

  # loop through all json cells
  for (i in 1:length(json$cells)){
    if (json$cells[[i]]$cell_type != "markdown"){ next }
    temp <- json$cells[[i]]$source
    temp2 <- replace.tags(temp,lab_text)
    json$cells[[i]]$source <- temp2
  }

  json2 <- RJSONIO::toJSON(json, asIs = F)

  outfile <- file.path(outdir,file)
  print(sprintf("Writing output to %s",outfile))  
  write(json2,outfile)

  cat("Parsed results in: ", outfile, "\n")
  cat("Please run convert_scanpy_labs.sh to execute notebooks and convert to html\n")
#  if (file == "scanpy_07_spatial.ipynb"){
#      print("scanpy_07_spatial was parsed, please run convert_scanpy_spatial.sh to execute notebook and convert to html")
#  }else{    
#    # render new notebook
#    render_nb <- sprintf("jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=1000 --inplace %s",outfile)
#    print("Convert to notebook...")
#    system(render_nb)
#
#    out_html <- sub(".ipynb",".html",outfile)  
#    render_html <- sprintf("jupyter nbconvert  --to html_toc --ExecutePreprocessor.timeout=1000  %s", outfile)
#    print("Convert to html...")
#    system(render_html)
#  }      
}
}

#jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=360 --inplace scanpy_01_qc.ipynb
#jupyter nbconvert --execute --to html_toc --ExecutePreprocessor.timeout=360 scanpy_01_qc.ipynb

######################################
# Scater / Seurat
######################################


pipelines <- c()
if (opt$type == "seurat"){
  pipelines <- c("seurat")
}else if (opt$type == "scater"){
  pipelines <- c("scater")
}else if (opt$type == "all"){
  pipelines <- c("seurat","scater")
}else if (opt$type != "scanpy"){
  stop(sprintf("Error! No pipeline %s, please specify one of scater, scanpy or seurat"))
}


for( pipeline in pipelines){
  outdir <- file.path(script_path,"compiled",pipeline)
  dir.create(outdir,recursive = T, showWarnings = F)
  
  indir <- pipeline
  if (opt$file == "all"){
    scripts <- dir(indir, pattern = "*.Rmd")
  }else {
    scripts <- opt$file
    if (!file.exists(file.path(script_path,indir,scripts))){
      stop(sprintf("File %s does not exist in dir %s",scripts,indir))
    }
  }
  
  for (script in scripts){
    infile <- file.path(script_path, pipeline, script)
    print(sprintf("Parsing data from %s",infile))
    lab <- readLines(infile)
    u <- grep("[#].*[_].*[:]",lab,value = T)[ grep("[#].*[_].*[:]",lab,value = T) 
                                              %in% sub(":.*",":",grep("[#].*[_].*[:]",lab_text,value = T)) ]
  
    t_lab <- lab
    for( i in u ){
      j <- grep(i,t_lab)
      temp2 <- sub(i,"",lab_text[grepl(i,lab_text)])
      t_lab <- c( t_lab[1:(j-1)] , temp2 , t_lab[(j+1):length(t_lab)] )
    }
  
    out_script <- file.path(outdir,script)
    print(sprintf("Writing output to %s",out_script))
    writeLines(t_lab, out_script)
    rmarkdown::render(out_script)
  }
}


