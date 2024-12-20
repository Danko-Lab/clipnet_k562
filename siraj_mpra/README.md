# QTL analysis of MPRA data from Siraj et al. 2024

This repository contains the code used to perform the QTL analysis of MPRA data from [Siraj et al. 2024](http://dx.doi.org/10.1101/2024.05.05.592437). The data used in this analysis was downloaded from the supplementary material of the paper.

Relevant addgene sequences:

https://www.addgene.org/109035/sequences/

## Process SI data

Use `process_SI.ipynb`.

## Generate reporter constructs

Run `create_reporter_constructs.py`. This script will generate the reporter constructs used in the MPRA experiment (may be worth confirming with the original authors), and random sequences to be use to potentially marginalize out minP activity.

```bash
python create_reporter_constructs.py media-3_oligos_snps_cleaned.tsv.gz ../data/mpra/k562_mpra_snps
python create_reporter_constructs.py media-3_oligos_snps_cleaned.tsv.gz ../data/mpra/k562_mpra_snps --dinuc_shuffles 1
python create_reporter_constructs.py media-3_oligos_snps_cleaned.tsv.gz ../data/mpra/k562_mpra_snps_2114 --procapnet
python create_reporter_constructs.py media-3_oligos_snps_cleaned.tsv.gz ../data/mpra/k562_mpra_snps_2114 --dinuc_shuffles 1 --procapnet
```

## Now predict on these sequences

```bash
conda activate clipnet

cd /home2/ayh8/clipnet/
python predict_ensemble.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref.h5 \
    --model_dir /home2/ayh8/clipnet_k562/models/clipnet_k562/ \
    --gpu 1
python predict_ensemble.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt.h5 \
    --model_dir /home2/ayh8/clipnet_k562/models/clipnet_k562/ \
    --gpu 0

for i in {1..9}; do
    python predict_individual_model.py \
        ../clipnet_k562/models/clipnet_k562/fold_${i}.h5 \
        ../clipnet_k562/data/mpra/k562_mpra_snps_ref.fa.gz \
        ../clipnet_k562/data/mpra/k562_mpra_snps_ref_fold_${i}.h5 \
        --gpu 0;
done

for i in {1..9}; do
    python predict_individual_model.py \
        ../clipnet_k562/models/clipnet_k562/fold_${i}.h5 \
        ../clipnet_k562/data/mpra/k562_mpra_snps_alt.fa.gz \
        ../clipnet_k562/data/mpra/k562_mpra_snps_alt_fold_${i}.h5 \
        --gpu 0;
done

python predict_ensemble.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_shuffle1.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_shuffle1.h5 \
    --model_dir /home2/ayh8/clipnet_k562/models/clipnet_k562/ \
    --gpu 1
```

## Predict using reference-trained CLIPNET

```bash
conda activate clipnet

cd ~/github/clipnet/
python predict_ensemble.py \
    ../data/k562_mpra/k562_mpra_snps_ref.fa.gz \
    ../predictions/k562/mpra/k562_mpra_snps_ref_ref_model.h5 \
    --model_dir ../clipnet_k562/models/clipnet_k562_reference/ \
    --gpu 0
python predict_ensemble.py \
    ../data/k562_mpra/k562_mpra_snps_alt.fa.gz \
    ../predictions/k562/mpra/k562_mpra_snps_alt_ref_model.h5 \
    --model_dir ../clipnet_k562/models/clipnet_k562_reference/ \
    --gpu 0

for i in {1..9}; do
    python predict_individual_model.py \
        ../clipnet_k562/models/clipnet_k562_reference/fold_${i}.h5 \
        ../data/k562_mpra/k562_mpra_snps_ref.fa.gz \
        ../predictions/k562/mpra/k562_mpra_snps_ref_fold_${i}_ref_model.h5 \
        --gpu 0;
done

for i in {1..9}; do
    python predict_individual_model.py \
        ../clipnet_k562/models/clipnet_k562_reference/fold_${i}.h5 \
        ../data/mpra/k562_mpra_snps_alt.fa.gz \
        ../predictions/k562/mpra/k562_mpra_snps_alt_fold_${i}_ref_model.h5 \
        --gpu 0;
done
```

## Predict using ProCapNet

```bash
conda activate procapnet

cd /home2/ayh8/clipnet_k562/siraj_mpra/
python predict_ensemble_procapnet.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_2114_ref.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_2114_ref.h5 \
    --model_dir /home2/ayh8/clipnet_k562/models/procapnet_k562/ \
    --gpu 0

python predict_ensemble_procapnet.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_2114_alt.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_2114_alt.h5 \
    --model_dir /home2/ayh8/clipnet_k562/models/procapnet_k562/ \
    --gpu 1
```

## Calculate predicted SNP effects

```bash
python calculate_mpra_predictions.py
```

## Create Enformer predictions

```bash
python predict_enformer_alt.py 
python predict_enformer_ref.py 
```

## Calculate Enformer SNP effects

```python
import h5py
import numpy as np
```

## Run DeepSHAP

```bash
conda activate clipnet

cd /home2/ayh8/clipnet/
allele=alt
mode=quantity
python calculate_deepshap.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_${allele}_fold0.fa \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_${allele}_fold0_deepshap_${mode}.npz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_${allele}_fold0_ohe.npz \
    --model_dir /home2/ayh8/clipnet_k562/models/clipnet_k562/ \
    --mode $mode \
    --gpu 0
```

## Run DeepSHAP for ProCapNet

```bash
conda activate procapnet

cd /home2/ayh8/clipnet_k562/siraj_mpra
allele=ref
mode=counts
for fold in {4..6}; do 
    python calculate_deepshap_procapnet.py \
        /home2/ayh8/clipnet_k562/models/procapnet_k562/fold_${fold}/ \
        /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_2114_${allele}_fold0.fa \
        /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_2114_${allele}_fold0_model${fold}_deepshap_${mode}.npz \
        --mode $mode \
        --gpu 0;
done
```
