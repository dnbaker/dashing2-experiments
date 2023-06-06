#!/bin/bash

# For a M1 Mac, where intel needs to be forced
docker build --tag dashing2-experiments --platform=linux/amd64 .

#docker run -it --platform=linux/amd64 fedora:39
