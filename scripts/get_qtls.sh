#!/bin/bash

zcat $1 | \
    awk -F, '$8=="TRUE"' | \
    sed -e '/chrX/d' -e '/chrY/d' | \
    awk -F, '{OFS="\t"}{print $1,$2,$2+1}'