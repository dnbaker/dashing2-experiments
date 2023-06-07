#!/bin/bash

set -ex

docker tag dashing2-experiments benlangmead/dashing2-experiments
docker push benlangmead/dashing2-experiments
