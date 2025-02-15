#!/bin/bash

n=$(nproc --all)

echo snakemake --printshellcmds --reason --jobs "$((n / 2))" --resources load=100
snakemake --printshellcmds --reason --jobs "$((n / 2))" --resources load=100