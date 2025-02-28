---
layout: default
title:  'FAQ'
---
<style>
h1, .h1, h2, .h2, h3, .h3, h4, .h4 { margin-top: 50px }
p.caption {font-size: 0.9em;font-style: italic;color: grey;margin-right: 10%;margin-left: 10%;text-align: justify}
</style>

### <img border="0" src="https://www.svgrepo.com/show/83019/faq-button.svg" width="40" height="40"> FAQ
***

<br/>

Below you can find some common error and problems you might face either during installations or during execution of softwares.


<details>
<summary><b>Command line developer tools not found (OSX)</b></summary>
<p>

If you don't yet have Mac OSX command line developer tools, please install it using:

```
xcode-select --install
```

</p>
</details>

<br/>


<details>
<summary><b>Error - umap-learn not found</b></summary>
<p>

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

</p>
</details>

<br/>







***
