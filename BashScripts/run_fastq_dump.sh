#!/usr/bin/env bash
for filename in /home/onur/workspace/rnaseq/data/*.sra; do
    /home/onur/workspace/rnaseq/student_tools/sratoolkit.current-ubuntu64/sratoolkit.2.11.2-ubuntu64/bin/fastq-dump $filename
done
