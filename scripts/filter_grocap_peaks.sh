#!/usr/bin/bash

cut -f 4- $1 | awk '$NF < 1000 {print}' | cut -f -3 | sort -k1,1 -k2,2n | uniq > $2
