#/usr/bin/env bash

awk '{sub(/[^,]*/,"");sub(/,/,"")} 1' $1