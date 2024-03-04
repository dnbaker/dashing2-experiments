#!/bin/bash

PLAT=aarch64

if [[ ! -d "dashing2-${PLAT}" ]] ; then
  git clone --recursive https://github.com/dnbaker/dashing2 -- "dashing2-${PLAT}"
fi
nm=$(grep docker build-${PLAT}.sh  | awk '{print $6}')

echo "When inside, run:"
echo "cd /code && make 2>&1 | tee log-${PLAT}.txt"

docker run -it -v`pwd`/dashing2-${PLAT}:/code --platform linux/${PLAT} ${nm}
