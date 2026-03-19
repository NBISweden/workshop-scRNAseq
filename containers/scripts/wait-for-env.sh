#!/bin/bash
# /usr/local/bin/wait-for-env

if [ ! -f /var/log/pixi/install.status ]; then
    echo "Waiting for pixi install to complete..."
    tail -f /var/log/pixi/install.log &
    TAIL_PID=$!
    while [ ! -f /var/log/pixi/install.status ]; do sleep 1; done
    kill $TAIL_PID
fi

if [ "$(cat /var/log/pixi/install.status)" = "ready" ]; then
    echo "Environment is ready!"
else
    echo "Environment install failed. Check /var/log/pixi/install.log"
    exit 1
fi
