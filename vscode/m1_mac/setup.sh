#!/bin/bash

set -ex

DR="${1}"
CONTAINING_DR=$(dirname "${DR}")
CLONE_DR=$(basename "${DR}")

pushd "${CONTAINING_DR}"
test ! -d "${CLONE_DR}"
git clone --recursive https://github.com/dnbaker/dashing2 -- "${CLONE_DR}"
popd

cp Dockerfile "${DR}"/
cp -r .devcontainer "${DR}"/
