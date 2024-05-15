#!/usr/bin/bash

awk '/>/{sub(">","&"FILENAME"_");sub(/\.fna/,x)}1' $1