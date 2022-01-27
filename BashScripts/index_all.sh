#!/usr/bin/env bash

for INFILE in "$@"
do
   samtools index $INFILE
done
