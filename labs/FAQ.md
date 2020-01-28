# Error - UMAP not found

If your R does not find the correct python version, it will complain that `umap-learn` is not installed and ask you to install it. Here are some tips on how to find the correct python version that was installed in the conda environment.

###  Try selecting the correct conda env in R:

In this example the conda env is named `sc_course`.

      library(reticulate)
      reticulate:use_conda("sc_course")

Then check what python you have in R:

     reticulate::py_config()
     # should read at top:
     python:         /Users/asbj/miniconda3/envs/sc_course/bin/python


If that still is not right, you may have an r-reticulate python installation as well and need to perform the steps below. 

### Restart R and select python version

OBS! Before doing anything else you need to select python version.

First, find out what path you have to your conda python (In terminal):

       which python
       /Users/asbj/miniconda3/envs/sc_course/bin/python

Then in R (after restarting):

     reticulate::use_python("/Users/asbj/miniconda3/envs/sc_course/bin/python", required=T)


Then check again with py_config if correct version of python is used:

     reticulate::py_config()

If you have the correct version now, you should be able to run UMAP without issues.


