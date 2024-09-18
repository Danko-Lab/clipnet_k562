# K562 PROcap data processing

This directory contains scripts for processing PRO-cap data from K562 cells into formats that can be used for training CLIPNET models.

## Download data

```bash
cd ./download_data
bash download.sh
```

## Get windows (including calling peaks)

```bash
cd ./get_window
bash make.sh
cd ../
```

## Process PRO-cap

```bash
cd ./procap
bash make.sh
cd ../
```

## Process sequence

```bash
cd ./sequence
bash make.sh
cd ../
```
