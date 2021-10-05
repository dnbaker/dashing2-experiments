#!/bin/bash
for f in `cat $1`;
do
    cat $f >> $2;
done
