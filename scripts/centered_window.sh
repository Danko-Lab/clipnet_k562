#!/bin/bash

zgrep -v 'chrom' $1 | grep -v 'NaN' | awk '{OFS="\t"}{print $1, $2-500,$2+500}'