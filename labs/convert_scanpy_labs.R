# Convert Scanpy labs to include tag text
# render as new ipynb and html.

# OBS! Requires jupyter_contrib_nbextensions to be installed.

# input is either --file all, of --file filename

library(optparse)
library(RJSONIO)

#parse input
option_list <- list( 
  make_option(c("-f", "--file"), type = "character", default="all",
              help="Notebook to convert, if or leave as all to run for all notebooks [all]")
)
opt <- parse_args(OptionParser(option_list=option_list))


#define path to this script 
initial.options <- commandArgs(trailingOnly = FALSE)
script_path <- dirname(sub("--file=","",initial.options[grep("--file=",initial.options)]))
setwd(script_path)

#opt <- list()
#opt$file <- "scanpy_04_clustering.ipynb"
#script_path <- "."

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

# read the ref text
lab_text <- readLines(paste0(script_path,"/knit_labs.Rmd"))

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

    # render new notebook
  render_nb <- sprintf("jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=360 --inplace %s",outfile)
  print("Convert to notebook...")
  system(render_nb)

  out_html <- sub(".ipynb",".html",outfile)  
  render_html <- sprintf("jupyter nbconvert --execute --to html_toc --ExecutePreprocessor.timeout=360  %s", outfile)
  print("Convert to html...")
  system(render_html)
}


#jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=360 --inplace scanpy_01_qc.ipynb
#jupyter nbconvert --execute --to html_toc --ExecutePreprocessor.timeout=360 scanpy_01_qc.ipynb





