#!/bin/bash

docker run -it --cap-add=SYS_PTRACE -v`pwd`:/code --platform=linux/amd64 dashing2-experiments
