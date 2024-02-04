#!/bin/sh

SIF=$1

export APPTAINERENV_USER=$(id -un)

## Get unused socket per https://unix.stackexchange.com/a/132524
## Tiny race condition between the python & singularity commands
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
cat 1>&2 <<END


        *************************************************
        *                                               *
        *  IMPORTANT: Do not close or exit this shell!  *
        *                                               *
        *************************************************

1. SSH tunnel from your workstation using the following command:

   ssh -N -L 8888:$(hostname):${PORT} ${APPTAINERENV_USER}@rackham.uppmax.uu.se

   point your web browser to http://localhost:8888/lab

2. Log in to Jupyter using the password:
   
   scrnaseq

When done using Jupyter, terminate the job by:

1. Shut down all kernels.
2. Issue the following command on the login node:

   CTRL-C (confirm 'y')

END

apptainer exec --bind ${PWD}:/run/user ${SIF} \
    /usr/local/bin/start.sh \
    jupyter lab \
        --no-browser \
        --ip=0.0.0.0 \
        --port=${PORT} \
        --ServerApp.password='argon2:$argon2id$v=19$m=10240,t=10,p=8$v6MTLSoeu/3mJPHOmiZ1sw$pr5VMmGV7zeOd2YZWu8lgP1lSBMtVeg/Mrj2XznRPEY' \
        --ServerApp.token='' \
        --allow-root

printf '\njupyter exited\n\n' 1>&2
