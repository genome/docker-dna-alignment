#!/bin/bash

set -o pipefail
set -o errexit

/usr/bin/java -Xmx4g -jar /opt/picard/picard.jar SamToFastq I=$1 INTERLEAVE=true INCLUDE_NON_PF_READS=true FASTQ=/dev/stdout | /usr/local/bin/bwa mem -K 100000000 -t $4 -Y -p -R "$2" $3 /dev/stdin | /usr/local/bin/samblaster -a --addMateTags | /opt/samtools/bin/samtools view -b -S /dev/stdin
