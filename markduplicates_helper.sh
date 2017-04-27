#!/bin/bash

set -o pipefail
set -o errexit

/usr/bin/java -Xmx16g -jar /opt/picard/picard.jar MarkDuplicates I=$1 O=/dev/stdout ASSUME_SORT_ORDER=queryname METRICS_FILE=mark_dups_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT | /usr/bin/sambamba sort -t $2 -m 18G -o $3 /dev/stdin
