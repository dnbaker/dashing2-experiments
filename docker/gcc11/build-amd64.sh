#!/bin/sh

TOOL=$(pwd | awk -v FS='/' '{print $NF}')
docker build --platform linux/amd64 -t dashing2-${TOOL}-amd64 .
