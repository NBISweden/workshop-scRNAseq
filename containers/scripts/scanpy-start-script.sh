#!/bin/bash
set -e

HOME="/home/jovyan"

pixi install --frozen --manifest-path ${HOME}/pixi.toml > /var/log/pixi/install.log 2>&1 && \
    echo "ready" > /var/log/pixi/install.status || \
    echo "failed" > /var/log/pixi/install.status &

function make_scanpy_kernel() (
    mkdir -p ${HOME}/.local/share/jupyter/kernels/scanpy
    cat <<EOF > ${HOME}/.local/share/jupyter/kernels/scanpy/kernel.json
{
  "argv": [
    "pixi", "run", "--frozen", "--manifest-path", "${HOME}/pixi.toml", "python",
    "-Xfrozen_modules=off",
    "-m",
    "ipykernel_launcher",
    "-f",
    "{connection_file}"
  ],
  "display_name": "scanpy",
  "language": "python",
  "metadata": {
    "debugger": true
  },
  "kernel_protocol_version": "5.5"
}
EOF
)

make_scanpy_kernel

exec jupyter lab \
    --no-browser \
    --ip=0.0.0.0 \
    --port=8888 \
    --PasswordIdentityProvider.hashed_password='argon2:$argon2id$v=19$m=10240,t=10,p=8$v6MTLSoeu/3mJPHOmiZ1sw$pr5VMmGV7zeOd2YZWu8lgP1lSBMtVeg/Mrj2XznRPEY' \
    --IdentityProvider.token='' \
    --MappingKernelManager.kernel_info_timeout=300
    --allow-root
