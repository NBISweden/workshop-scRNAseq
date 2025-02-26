#!/bin/bash

source /opt/conda/bin/activate scanpy

exec jupyter lab \
    --no-browser \
    --ip=0.0.0.0 \
    --port=8888 \
    --PasswordIdentityProvider.hashed_password='argon2:$argon2id$v=19$m=10240,t=10,p=8$v6MTLSoeu/3mJPHOmiZ1sw$pr5VMmGV7zeOd2YZWu8lgP1lSBMtVeg/Mrj2XznRPEY' \
    --IdentityProvider.token='' \
    --allow-root
