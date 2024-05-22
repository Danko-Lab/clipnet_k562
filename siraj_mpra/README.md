# QTL analysis of MPRA data from Siraj et al. 2024

This repository contains the code used to perform the QTL analysis of MPRA data from [Siraj et al. 2024](http://dx.doi.org/10.1101/2024.05.05.592437). The data used in this analysis was downloaded from the supplementary material of the paper.

Relevant addgene sequences:

https://www.addgene.org/109035/sequences/

## Process SI data

Use `process_SI.ipynb`.

## Generate reporter constructs

Run `create_reporter_constructs.py`. This script will generate the reporter constructs used in the MPRA experiment (may be worth confirming with the original authors), and random sequences to be use to potentially marginalize out minP activity.

```bash
python create_reporter_constructs.py media-3_oligos_snps.tsv.gz ../data/mpra/k562_mpra_snps
python create_reporter_constructs.py media-3_oligos_snps.tsv.gz ../data/mpra/k562_mpra_snps --dinuc_shuffles 1
python create_reporter_constructs.py media-3_oligos_snps.tsv.gz ../data/mpra/k562_mpra_snps_2114 --procapnet
```

## Now predict on these sequences

```bash
conda activate clipnet

cd /home2/ayh8/clipnet/
python predict_ensemble.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref.h5 \
    --model_dir /home2/ayh8/clipnet_k562/models/clipnet_k562/ \
    --gpu --use_specific_gpu 1
python predict_ensemble.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt.h5 \
    --model_dir /home2/ayh8/clipnet_k562/models/clipnet_k562/ \
    --gpu --use_specific_gpu 0
python predict_ensemble.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_shuffle1.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_shuffle1.h5 \
    --model_dir /home2/ayh8/clipnet_k562/models/clipnet_k562/ \
    --gpu --use_specific_gpu 0

cd ../

```
