#!/bin/bash

grep -v 'chrom' $1 | \
    grep -v 'NaN' | \
    awk '{OFS="\t";OFMT="%f"}{print $1, int(($2+$3)/2)-500, int(($2+$3)/2)+500}'