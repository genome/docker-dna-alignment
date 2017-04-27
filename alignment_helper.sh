#!/bin/bash

set -o pipefail
set -o errexit

/usr/local/bin/bwa mem -K 100000000 -t $5 -Y -R "$1" $2 $3 $4 | /usr/local/bin/samblaster -a --addMateTags | /opt/samtools/bin/samtools view -b -S /dev/stdin
