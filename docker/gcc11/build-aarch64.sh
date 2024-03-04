#!/bin/sh

TOOL=$(pwd | awk -v FS='/' '{print $NF}')
docker build --platform linux/aarch64 -t dashing2-${TOOL}-aarch64 .
