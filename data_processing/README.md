# K562 PROcap data processing

This directory contains scripts for processing PRO-cap data from K562 cells into formats that can be used for training CLIPNET models. The data can be downloaded and processed using the following pipelines. Following the execution of these pipelines, the data can then be used for fine tuning and evaluating CLIPNET models. NOTE: These scripts contain hardcoded path variables at the top of each `Snakefile`. To use them, please modify these variables so they point to appropriate paths on your system.

## Download raw data from ENCODE

The K562 PRO-cap data and ENCODE reference genome copy can be downloaded using the bash script in `download_data`:

```bash
cd ./download_data
bash download.sh
```

## Get windows (including calling peaks)

The snakemake pipeline in `get_window` will call peaks & generate 1kb windows for training.

```bash
cd ./get_window
bash make.sh
```

## Process PRO-cap

The snakemake pipeline in `procap` will RPM normalize the PRO-cap data and extract the peak regions for training.

```bash
cd ./procap
bash make.sh
```

## Process sequence

The snakemake pipeline in `sequence` will extract the sequence around peak regions.

```bash
cd ./sequence
bash make.sh
```

## PRO-seq data processing

All other pipelines are for processing data for K562 PRO-seq and pausing models that are still in-progress (if you have a need for such a model or want to explore sequence -> pausing relationships, please email me to discuss).
