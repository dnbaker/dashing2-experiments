#!/bin/bash

set -ex

find . -name '*.sig'
find . -name '*.opss' | xargs rm -f
find . -name '*.ss' | xargs rm -f
find . -name 'fastani.tsv' | xargs rm -f
