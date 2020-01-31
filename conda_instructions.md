---
layout: default
title:  'Schedule - scRNAseq course'
---

#### <img border="0" src="https://hackernoon.com/hn-images/1*rW03Wtue71AKfxnx6XN_iQ.png" width="50" height="50"> Conda Instructions
***

In this workshop you will use conda environments to run the exercises. This is because conda environments allow all students to have the save computing environment, i.e. package versions. This enforces reproducibility for you to run this material without the need to re-install or change your local versions. See and graphical example below:


<img border="0" src="https://nbisweden.github.io/excelerate-scRNAseq/logos/conda_illustration.png" width="400">


[Conda environments](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) are a self-contained directory that
you can use in order to reproduce all your results. Two of the required software are not available as Conda packages, please see the separate [instructions for installing SingleR and CHETAH](https://raw.githubusercontent.com/NBISweden/excelerate-scRNAseq/master/notes_installation.txt).

Briefly, you need to:  

1. Install Conda and download the `.yml` file
2. Create and activate the environment
3. Deactivate the environment after running your analyses

You can [read more](https://nbis-reproducible-research.readthedocs.io/en/latest/conda/) about Conda environments and other important concepts to help you make your research reproducible.

<br/>

<br/>

##### <img border="0" src="https://www.svgrepo.com/show/4795/installation-symbol.svg" width="20" height="20"> Install Conda and download the environment file
***

You should start by installing Conda. We suggest installing either Miniconda3 (NOT Anaconda). After [installing Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), download the course [Conda file](labs/environment_r.yml) and put it in your working folder.


###### **On MacOSX**

```
curl -o Miniconda3-latest-MacOSX-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
sh Miniconda3-latest-MacOSX-x86_64.sh
```

Follow the instructions on screen replying `yes` when necessary. Restart your terminal window to apply modifications.


###### **On Ubuntu**

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
```

Follow the instructions on screen replying `yes` when necessary. Restart your terminal window to apply modifications.


###### **On Windows10**

Several packages are not available for Windows. However, on windows10 we can run a Ubuntu subsystem to overcome this issue. Please follow the instructions `Alternative option on Windows (WLS)` below to install it.

<br/>

<br/>

##### <img border="0" src="https://www.svgrepo.com/show/4795/installation-symbol.svg" width="20" height="20"> Create a Conda environment from file
***

To download the `environment_r.yml` file using the command on Terminal:

```
#Ubuntu
wget https://nbisweden.github.io/workshop-scRNAseq/labs/environment_r.yml

#MacOSX
curl -o environment_r.yml https://nbisweden.github.io/workshop-scRNAseq/labs/environment_r.yml
```

After this, you should have a file named `environment_r.yml` in your directory (it does not matter where, you can save on Downloads folder for example). Next, type:

```
conda env create -p scRNAseq2020 -f environment_r.yml
```

Several messages will show up on your screen and will tell you about the installation process. This may take a few minutes depending on how many packages are to be installed.

```
##Collecting package metadata: done
##Solving environment: done
##
##Downloading and Extracting Packages
##libcblas-3.8.0       | 6 KB      | ############################################################################# | 100%
##liblapack-3.8.0      | 6 KB      | ############################################################################# | 100%
##...
##Preparing transaction: done
##Verifying transaction: done
##Executing transaction: done
```

<br/>

<br/>

##### <img border="0" src="https://www.svgrepo.com/show/4795/installation-symbol.svg" width="20" height="20"> Activate the environment
***

Once the environment is created, we need to activate it in order to use the softwares and packages inside it. To activate an environment type:

```
source activate scRNAseq2020
```

From this point on you can run any of the contents from the course. For instance, you can directly launch RStudio by
typing `rstudio`. Here it is important to add the `&` symbol in the end to be able to use the command line at the same time if needed. You can open other files from Rstudio later as well.

```
rstudio ./labs/compiled/my_script.Rmd &
```

Similarly, you can open python notebooks by typing:

```
jupyter notebook ./labs/scapy/01_qc.ipynb &
```

<br/>

<br/>

##### <img border="0" src="https://www.svgrepo.com/show/4795/installation-symbol.svg" width="20" height="20"> Deactivate the environment
***

After you've ran all your analyses, you can deactivate the environment by typing:

```
conda deactivate
```

<br/>

<br/>



##### <img border="0" src="https://upload.wikimedia.org/wikipedia/commons/5/5f/Windows_logo_-_2012.svg" width="20" height="20"> Alternative option on Windows (WLS)
***

Unfortunately, not all packages available on conda are compatible with windows machines. The good news is that is changed on windows10, in which they offer native linux support via the Windows Subsystem for Linux (WSL2). This allows you to run linux/bash commands from within windows without the need of a virtual machine nor a dual-boot setup (i.e. having 2 operational systems). However, WSL does not offer a complete support for graphical interfaces (such as RStudio in our case), so we need an additional steps to make that happen.

1. On Windows10, install the WSL if you don't have it. Follow the instructions here:
[https://docs.microsoft.com/en-us/windows/wsl/install-win10](https://docs.microsoft.com/en-us/windows/wsl/install-win10)

2. Once you have that installed, you can download and install MobaXterm (which is the enhanced terminal with graphical capacity):
[https://mobaxterm.mobatek.net](https://mobaxterm.mobatek.net)

3. Inside MobaXterm, you will probably will see that your WSL is already listed on the left panel as an available connection. Just double-click it and you will be accessing it via MobaXterm. If by any chance you don't see it there, close MobaXterm and go to the WSL terminal, because probably the WSL is not allowing SSH connections. You can follow this [link](https://www.illuminiastudios.com/dev-diaries/ssh-on-windows-subsystem-for-linux/) for the instructions on how to do it. You need to complete until the step `Start or restart the SSH service`, while the further steps are optional, but might be useful.

4. Inside MobaXterm, download Conda with the command:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

5. Inside MobaXterm, type the commands below to install Conda. Follow the instructions for the installation there.
```
cd ~/Downloads
sh Miniconda3-latest-Linux-x86_64.sh
```

6. Inside MobaXterm, type the commands below to install the X-server graphical packages that will be used to launch RStudio.
[https://docs.anaconda.com/anaconda/install/linux/](https://docs.anaconda.com/anaconda/install/linux/)
```
sudo apt-get update
sudo apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
```

7. Close Terminal to apply the CONDA updates. Then you can create a course folder, download the environment file and create the environment:
```
mkdir ~/Desktop/course
cd ~/Desktop/course
wget https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/master/labs/environment_r.yml
conda env create -n scRNAseq2020 -f environment_r.yml
```

8. You can then follow the instructions above to activate/deactivate the environment.
```
conda activate scRNAseq2020
rstudio &
```

<br/>

<br/>


##### <img border="0" src="https://upload.wikimedia.org/wikipedia/commons/d/d5/Virtualbox_logo.png" width="20" height="20"> Alternative option (VIRTUALBOX)
***

If by any means you see that the installations are not working as it should on your computer, you can try to create a virtual machine to run UBUNTU and install everything there.

1. Download and install on your machine VIRTUALBOX
[https://www.virtualbox.org](https://www.virtualbox.org)

2. Download the ISO disk of UBUNTU
[https://ubuntu.com/download/desktop](https://ubuntu.com/download/desktop)

3. On VIRTUALBOX, click on `Settings` (yellow engine) > `General` > `Advanced` and make sure that both settings **Shared Clipboard** and **Drag'n'Drop** are set to `Bidirectional`.

4. Completely close VIRTUALBOX and start it again to apply changes.

5. On VIRTUALBOX, create a machine called Ubuntu and add the image above
- set the memory to the maximum allowed in the GREEN bar
- set the hard disk to be dynamic allocated
- all other things can be default

6. Proceed with the Ubuntu installation as recommended. You can set to do "Minimal Installation" and deactivate to get updates during installation.


7. Inside Ubuntu, open TERMINAL and type the commands below to install the X-server graphical packages that will be used to launch RStudio.
[https://docs.anaconda.com/anaconda/install/linux/](https://docs.anaconda.com/anaconda/install/linux/)
```
sudo apt-get update
sudo apt-get install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6
```

8. Inside UBUNTU, Download conda:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

9. Inside UBUNTU, open the TERMINAL and type the commands below. Follow the instructions for the installation there.
```
cd ~/Downloads
sh Miniconda3-latest-Linux-x86_64.sh
```

10. Close Terminal to apply the CONDA updates. Then you can create a course folder, download the environment file and create the environment:
```
mkdir ~/Desktop/course
cd ~/Desktop/course
wget https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/master/labs/environment_r.yml
conda env create -n scRNAseq2020 -f environment_r.yml
```

11. You can then follow the instructions above to activate/deactivate the environment.
```
conda activate scRNAseq2020
rstudio &
```
