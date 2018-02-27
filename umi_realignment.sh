#!/bin/bash

set -o pipefail
set -o errexit
set -o nounset

if [ $# -ne 3 ]
then
    echo "Usage: $0 unaligned.bam reference.fa cores"
    exit 1
fi

BAM="$1"
REFERENCE="$2"
CORES="$3"

/usr/bin/java -Xmx1g -jar /opt/picard/picard.jar SamToFastq \
    "I=$BAM" INCLUDE_NON_PF_READS=true F=/dev/stdout \
    INTERLEAVE=true \
 | /usr/local/bin/bwa mem -p -t "$CORES" "$REFERENCE" /dev/stdin \
 | /usr/bin/java -Xmx4g -jar /opt/picard/picard.jar MergeBamAlignment \
    "UNMAPPED=$BAM" ALIGNED=/dev/stdin OUTPUT=realigned.bam "REFERENCE_SEQUENCE=$REFERENCE" \
    CLIP_ADAPTERS=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true \
    EXPECTED_ORIENTATIONS=FR MAX_GAPS=-1 SORT_ORDER=coordinate \
    ALIGNER_PROPER_PAIR_FLAGS=false \
    ATTRIBUTES_TO_RETAIN=X0 ATTRIBUTES_TO_RETAIN=ZS ATTRIBUTES_TO_RETAIN=ZI \
    ATTRIBUTES_TO_RETAIN=ZM ATTRIBUTES_TO_RETAIN=ZC ATTRIBUTES_TO_RETAIN=ZN \
    ATTRIBUTES_TO_REVERSE=ad ATTRIBUTES_TO_REVERSE=bd ATTRIBUTES_TO_REVERSE=cd \
    ATTRIBUTES_TO_REVERSE=ae ATTRIBUTES_TO_REVERSE=be ATTRIBUTES_TO_REVERSE=ce
