#!/bin/sh

docker buildx build --platform linux/amd64 -t dashing2-gcc9 .
