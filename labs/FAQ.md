# Frequently Asked Questions

Here we will try to fill up questions that may arise during the exercises that we feel would be useful for many course participants to have a look at.

## maximal number of DLLs reached...

Many of the scRNAseq packages have a heap of different dependencies. There is a limit to how many libraries you can load at the same time, and when loading a new library you may run into an error message like:

     	Error: package or namespace load failed for ‘Seurat’ in dyn.load(file, DLLpath = DLLpath, ...):
      	unable to load shared object '/Library/Frameworks/R.framework/Versions/3.4/Resources/library/tclust/libs/tclust.so':
      `	maximal number of DLLs reached...
      Error: package or namespace load failed for ‘Seurat’ in dyn.load(file, DLLpath = DLLpath, ...):

If this happens the easiest is to start a new R session and start again from scratch. Or to remove most of your loaded packages. You can also modify the maximal number of DLLs with the environmental variable R_MAX_NUM_DLLS before starting R to permit more loaded packages. 

On Mac  mofify the file /Library/Frameworks/R.framework/Resources/etc/Renviron and add R_MAX_NUM_DLLS=500 in the end. For other OS check https://stat.ethz.ch/R-manual/R-patched/library/base/html/Startup.html.

