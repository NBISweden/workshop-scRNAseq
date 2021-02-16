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
  <summary markdown="span">`Error - umap-learn not found`</summary>

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
  curl -o MacOSX10.9.sdk.tar.gz "https://github-production-release-asset-2e65be.s3.amazonaws.com/13597203/f0123b00-34ab-11ea-84b1-27ccc324f983?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAIWNJYAX4CSVEH53A%2F20210122%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20210122T065821Z&X-Amz-Expires=300&X-Amz-Signature=e11864967b0c9a1e1ba1121dbadd35bb3129ae1cd87bdf07b1a9965c731ae129&X-Amz-SignedHeaders=host&actor_id=22674952&key_id=0&repo_id=13597203&response-content-disposition=attachment%3B%20filename%3DMacOSX10.9.sdk.tar.xz&response-content-type=application%2Foctet-stream"

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


***
