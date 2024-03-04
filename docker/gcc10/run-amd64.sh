#!/bin/bash

set -ex

TOOL=$(pwd | awk -v FS='/' '{print $NF}')
PLAT=amd64

if [[ ! -d "dashing2-${PLAT}" ]] ; then
  git clone --recursive https://github.com/dnbaker/dashing2 -- "dashing2-${PLAT}"
fi

docker run -it -v`pwd`/dashing2-${PLAT}:/code --platform linux/${PLAT} \
  dashing2-${TOOL}-${PLAT} \
  /bin/bash -c "cd /code && make 2>&1 | tee log-${PLAT}.txt"

cp dashing2-${PLAT}/log-${PLAT}.txt .
