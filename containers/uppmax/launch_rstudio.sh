#!/bin/sh

SIF=$1

## Create temporary directory to be populated with directories to bind-mount in the container
## where writable file systems are necessary. Adjust path as appropriate for your computing environment.
workdir=$(python -c 'import tempfile; print(tempfile.mkdtemp())')

echo "container workdir: ${workdir}"

mkdir -p -m 700 ${workdir}/run ${workdir}/tmp ${workdir}/var/lib/rstudio-server
cat > ${workdir}/database.conf <<END
provider=sqlite
directory=/var/lib/rstudio-server
END

## Set OMP_NUM_THREADS to prevent OpenBLAS (and any other OpenMP-enhanced
## libraries used by R) from spawning more threads than the number of processors
## allocated to the job.
##
## Set R_LIBS_USER to a path specific to rocker/rstudio to avoid conflicts with
## personal libraries from any R installation in the host environment

cat > ${workdir}/rsession.sh <<END
#!/bin/sh
export OMP_NUM_THREADS=$(nproc --all)
export R_LIBS_USER=/usr/local/lib/R/site-library
exec /usr/lib/rstudio-server/bin/rsession "\${@}"
END

chmod +x ${workdir}/rsession.sh

export APPTAINER_BIND="${workdir}/run:/run,${workdir}/tmp:/tmp,${workdir}/database.conf:/etc/rstudio/database.conf,${workdir}/rsession.sh:/etc/rstudio/rsession.sh,${workdir}/var/lib/rstudio-server:/var/lib/rstudio-server"

## Do not suspend idle sessions.
## Alternative to setting session-timeout-minutes=0 in /etc/rstudio/rsession.conf
## https://github.com/rstudio/rstudio/blob/v1.4.1106/src/cpp/server/ServerSessionManager.cpp#L126
export APPTAINERENV_RSTUDIO_SESSION_TIMEOUT=0

export APPTAINERENV_USER=$(id -un)
export APPTAINERENV_PASSWORD="scrnaseq"
## Uncomment the next line for a randomly generated password
# export APPTAINERENV_PASSWORD=$(openssl rand -base64 15)

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

   ssh -N -L 8787:${HOSTNAME}:${PORT} ${APPTAINERENV_USER}@${HOSTNAME}

   and point your web browser to http://localhost:8787

2. log in to RStudio Server using the following credentials:

   user: ${APPTAINERENV_USER}
   password: ${APPTAINERENV_PASSWORD}

When done using RStudio Server, terminate the job by:

1. Exit the RStudio Session ("power" button in the top right corner of the RStudio window)
2. Issue the following command on the login node:

      CTRL-C

END

apptainer exec --cleanenv ${SIF} \
    /usr/lib/rstudio-server/bin/rserver --server-user ${APPTAINERENV_USER} --www-port ${PORT} \
            --auth-none=0 \
            --auth-pam-helper-path=pam-helper \
            --auth-stay-signed-in-days=30 \
            --auth-timeout-minutes=0 \
            --rsession-path=/etc/rstudio/rsession.sh

printf 'rserver exited' 1>&2
