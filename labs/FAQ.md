---
layout: default
title:  'FAQ'
---

{::options parse_block_html="true" /}

<style>
h1, .h1, h2, .h2, h3, .h3, h4, .h4 { margin-top: 50px }
p.caption {font-size: 0.9em;font-style: italic;color: grey;margin-right: 10%;margin-left: 10%;text-align: justify}
</style>

### <img border="0" src="https://www.svgrepo.com/show/83019/faq-button.svg" width="40" height="40"> FAQ
***

<br/>

Below you can find some common error and problems you might face either during installations or during execution of softwares.


<details>
  <summary markdown="span">`Command line developer tools not found` (OSX)</summary>

  If you don't yet have Mac OSX command line developer tools, please install it using:

  ```
  xcode-select --install
  ```

</details>


<details>
  <summary markdown="span">`Error - umap-learn not found, or other python packagees`</summary>

  If your R does not find the correct python version, it will complain that `umap-learn` is not installed and ask you to install it. Here are some tips on how to find the correct python version that was installed in the conda environment.

  <br/>

  **Try selecting the correct conda env in R**

  In this example the conda environment is named `scRNAseq2021`.
  ```
  library(reticulate)
  reticulate::use_conda("scRNAseq2021")
  ```

  Then check what python you have in R:
  ```
  reticulate::py_config()
  # should read at top:
  python:         /Users/asbj/miniconda3/envs/scRNAseq2021/bin/python
  ```

  If that still is not right, you may have an `r-reticulate` python installation as well and need to perform the steps below.

  <br/>

  **Restart R and select python version**

  OBS! Before doing anything else you need to select python version.

  First, find out what path you have to your conda python (in TERMINAL):
  ```
  which python
  /Users/asbj/miniconda3/envs/scRNAseq2021/bin/python
  ```

  Then in R (after restarting):
  ```
  reticulate::use_python("/Users/asbj/miniconda3/envs/scRNAseq2021/bin/python", required=T)
  ```

  Then check again with `py_config` if correct version of python is used:
  ```
  reticulate::py_config()
  ```

  If you have the correct version now, you should be able to run UMAP without issues.

</details>

<details>
  <summary markdown="span">`Unable to load stringi.so` (UNIX/Windows)</summary>

  You can install stringi in R using:

  ```
  install.packages('stringi')
  ```

</details>


<details>
  <summary markdown="span">`ERROR: Failed building wheel for gevent` / `MacOSX10.9.sdk missing` (MacOSX)</summary>

  This is a problem with the MacOSX compiler, in which conda is unable to find it.

  ```
  #Download MacOSX10.9.sdk from Github
  curl -o MacOSX10.9.sdk.tar.gz "https://github.com/phracker/MacOSX-SDKs/releases/download/11.3/MacOSX10.9.sdk.tar.xz"

  #extract
  sudo tar -xzf MacOSX10.9.sdk.tar.xz

  #copy
  sudo cp -r MacOSX10.9.sdk /opt/

  #give executable permissions
  sudo chmod -R a+rX /opt

  #Link the path where conda looks to where the file is
  ln -s /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk /opt/MacOSX10.9.sdk
  ```

</details>

<details>
  <summary markdown="span">`ERROR: option error has NULL value`</summary>

  This error happens when running code inline.
  
  One possible solution is to restart Rstudio and type.
  
  ```
  if(interactive()) { options(error = utils::recover)}
  ```
  
  Please try other solutions listed here: https://github.com/rstudio/rstudio/issues/4723
  
  If none of those work, you can click on the wheel engine symbol and check `Chunk output in console` 

</details>


<details>
  <summary markdown="span">`R crashes due to memory issues`</summary>

  If R crashes due to memory issues, it may be a good idea to increase the vector size `R_MAX_VSIZE`. Put in the file `.Renviron` either in your home directory or the folder you are launching Rstudio from:

  ```
  R_MAX_VSIZE=70Gb
  ```
  Or to whatever value matches your computer, the default size is 16Gb.


</details>



***
