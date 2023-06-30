#!/bin/bash

set -ex

./00wide.R
./01pivoted.R
./02augmented.R
./03filtered.R
./04sse.R
./06augmented_subset.R
./07filtered_subset.R
./20sedf.R
