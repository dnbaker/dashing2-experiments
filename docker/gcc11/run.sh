#!/bin/bash

if [[ ! -d dashing ]] ; then
  git clone --recursive https://github.com/dnbaker/dashing2 
fi
nm=$(grep docker build.sh  | awk '{print $4}')

echo "When inside, run:"
echo "cd /code && make -j4"

docker run -it -v`pwd`/dashing2:/code ${nm}
