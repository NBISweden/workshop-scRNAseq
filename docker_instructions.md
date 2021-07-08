# <img border="0" src="https://cdn1.iconfinder.com/data/icons/social-media-2106/24/social_media_social_media_logo_docker-512.png" width="40" height="40"> Docker

***

<br/>

In the last case scenario, if you are having problems installing R packages, please follow these instructions to setup a docker container with R and Rstudio. Then you will need to install packages as usual. See the [pre-course material](/single-cell_sib_scilifelab/precourse.md) for package installations.

<img border="0" src="https://quppler.com/wp-content/uploads/2019/03/DockerComponents.png" width="500">

<br/>

## <img border="0" src="https://www.svgrepo.com/show/4795/installation-symbol.svg" width="40" height="40"> Installing Docker

***

1.Install docker. You need permission to both install and use `sudo` commands.

* Mac: https://docs.docker.com/docker-for-mac/install/
* Windows: https://docs.docker.com/docker-for-windows/install/
* Ubuntu: https://docs.docker.com/install/linux/docker-ce/ubuntu/

2.Launch docker and got to `Preferences > Resources > Advanced`  and set `Memory=8gb, cpu=4, swp=1`. These will set the maximum values possible for your machine.

3.Now go to `Preferences > Security` and tick on "Allow connections from network clients", so that you can communicate with DockerServer.

4.Open a terminal shell.

5.You will need to use `docker-machine`. If your docker version does not automatically install `docker-machine` in terminal, please [follow these instructions](https://github.com/docker/machine/releases/tag/v0.16.2) to install it. Then create a Virtual Machine (VM) named "default" and set the amount of CPU and RAM available for you:

```bash
docker-machine create default
VBoxManage modifyvm default --cpus 4
VBoxManage modifyvm default --memory 8192
```

6.You can now start the machine, the environment (which will give you an IP address to connect to it!) and configure your shell to work with docker.

```bash
docker-machine start
docker-machine env
eval "$(docker-machine env default)"
```

If at any point from now on you get an error asking you to stop the container, you can do so with:

```bash
docker-machine stop
```

<br/>

## <img border="0" src="https://www.svgrepo.com/show/264/user.svg" width="40" height="40"> Using the SingleCellSchool container

***

Go to a folder where you can safelly download the course materials (to make them visible to the container). you can simply use 

```bash
git clone https://github.com/NBISweden/single-cell_sib_scilifelab
```

or

```bash
cd ~/
svn export -O https://github.com/NBISweden/single-cell_sib_scilifelab/trunk
mv trunk SingleCellSchool
cd SingleCellSchool
```

Now that you are setup with docker you can launch/fire the container with all packages from the course from a repository. If you already have a copy on your computer it will launch it, otherwise it will download it from the repositiory first and then launch it. `czarnewski/single_cell_school` is where the docker image is located.

```bash
docker run -d --rm -v $(pwd):/home/rstudio --memory=8g -p 8787:8787 -e PASSWORD=test czarnewski/single_cell_school
```
If no errors are thrown, you will be able to connect to your Rstudio machine using a webbrowser, just type the IP.

* http://192.168.99.100:8787 (this might change, check the output from the command `docker-machine env` on the step above)
* LOGIN=`rstudio`
* PASS=`test`

You can alternativelly check your docker container IP and replace the (192.168.99.100) above:

```bash
docker-machine ip
```

Your RStudio session should look like this now:

<img border="0" src="https://github.com/NBISweden/single-cell_sib_scilifelab/raw/master/logos/SingleCellSchool_Docker.png" width="500">

You should have for example the package `Seurat` installed and also be able to see the files from the course on your home folder. Just by navigating thought there you will be able to open the exercise files (.Rmd extension). 

You can save your container back to an image at any time by typing into the terminal. To run it again just repeat the `docker run` step above .

```bash
docker commit <container_ID> YOURusername/NEWimagename
```

You can stop a container from running with:

```bash
docker container stop <container_ID>
```


<br/>

## <img border="0" src="https://www.svgrepo.com/show/2270/pin-tool.svg" width="40" height="40"> Building a Docker Container from scratch

***

7.Now we can download a pre-made docker container containing both the latest version R and RStudio (rocker). Here we need to set a new password for the container, which we can set as "test". The `-p 8787:8787` indicates which ports are visible between your computer and the container for visualizing Rstudio. `--rm` will remove the container after use. `-e` sets the password.

```bash
docker run -d --rm --memory=8g -p 8787:8787 -e PASSWORD=test --name rstudio rocker/verse
```

If no errors are thrown, you will be able to connect to your Rstudio machine using a webbrowser, just type the IP.

* http://192.168.99.100:8787 (this might change, check the output from the command `docker-machine env` on the step above)
* LOGIN=`rstudio`
* PASS=`test`

You can alternativelly check your IP using:

```bash
docker-machine ip
```

8.You can now proceed with using Rstudio and installing packages as usual.

See the [pre-course material](/single-cell_sib_scilifelab/precourse.md)

9.To save your container status, you can commit the changes back to a new image using. If you plan on sharing your image, `YOURusername/NEWimagename` must be a valid Repository on your account on DockerHub (see below).

```bash
#Login into your DockerHub account
docker login --username=YOURusername --email=youremail@company.com

#Get container ID
docker container ls

#Create a commit to save alterations to a new image
docker commit <container_ID> YOURusername/NEWimagename
```

Now you can access your image container at any time by running.

```bash
docker run -d -rm -v $(pwd):/home/rstudio --memory=8g -p 8787:8787 -e PASSWORD=test czarnewski/single_cell_school
```

If no errors are thrown, you will be able to connect to your Rstudio machine using a webbrowser, just type the IP.

* http://192.168.99.100:8787 (this might change, check the output from the command `docker-machine env` on the step above)
* LOGIN=`rstudio`
* PASS=`test`

You can alternativelly check your docker container IP and replace the (192.168.99.100) above:

```bash
docker-machine ip
```

<br/>

## <img border="0" src="https://www.svgrepo.com/show/23574/share-with-cloud.svg" width="40" height="40"> Sharing your Docker Image

***

1. Create a account/repository in your DockerHub.

* Log in on https://hub.docker.com/
* Click on `+ Create Repository`.
* Choose a name (e.g. NEWimagename) and a description for your repository and click Create.
* Your image name from now on should be called `YOURusername/NEWimagename` (continue below to see how). 

You can apply a tag as of a version of your container. And then push it to your repository

```bash
docker tag <container_ID> YOURusername/NEWimagename:latest
docker push yourhubusername/gapminder_my_analysis
```

<br/>


## <img border="0" src="https://cdn4.iconfinder.com/data/icons/proglyphs-computers-and-development/512/Terminal-512.png" width="40" height="40"> Useful Commands and links

***

Some of these commands can be used in

```bash
#get container ip address
docker-machine ip

#list or remove docker containers
docker container ls
docker rm <container_id>

#list or remove docker images
docker image ls
docker rmi <image_id>

#runs a interactive (-ti) shell inside the container.
#"-v" mounts the current directory "$(pwd)" to the folder "/rstudio".
#You can use them separately
docker run -it -v $(pwd):/home/rstudio <docker_name>
```

Usefull Links:

* [https://ropenscilabs.github.io/r-docker-tutorial/](https://ropenscilabs.github.io/r-docker-tutorial/)
* [https://rmarkdown.rstudio.com/articles_integration.html](https://rmarkdown.rstudio.com/articles_integration.html)
* [https://github.com/rocker-org/rocker/wiki/How-to-save-data](https://github.com/rocker-org/rocker/wiki/How-to-save-data)

<br/>

## [Back to main](README.md)
