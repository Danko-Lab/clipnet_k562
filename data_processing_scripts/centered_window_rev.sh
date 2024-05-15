#!/bin/bash

awk '{OFS="\t"}{print $1, $3-500,$3+500}' $1