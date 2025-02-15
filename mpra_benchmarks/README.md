# Deep learning analysis of MPRA data from Siraj et al. 2024

This repository contains the code used to benchmark deep learning models on MPRA data from [Siraj et al. 2024](http://dx.doi.org/10.1101/2024.05.05.592437). The data used in this analysis was downloaded from the supplementary material of the paper.

## Process SI data

The SI from the paper are provided as xlsx files, which can be slightly annoying to parse. The notebook `process_SI.ipynb` contains code that downloads and transforms these xlsx files into TSVs for easier parsing. The two main TSVs I used are `media-3_oligos_snps_cleaned.tsv.gz` and `media-4-K562_allelic_mpra.tsv.gz`. These have been archived on [Zenodo]().

## Generate reporter constructs

Run `create_reporter_constructs.py`. This script will generate the reporter constructs used in the MPRA experiment (backbone sequences taken from [Addgene](https://www.addgene.org/109035/sequences/)). This script also has the option to insert dinucleotide shuffles of the enhancer constructs instead, which might be useful if you wanted to, say, marginalize out minP activity (i.e., what is the predicted regulatory activity of the minP/backbone with an inert insert), but this wasn't performed as part of our study. Using the `--procapnet` flag will extend the sequences out to 2114 bp (the input length for ProCapNet/BPNet) by including more of the reporter backbone.

```bash
python create_reporter_constructs.py media-3_oligos_snps_cleaned.tsv.gz ../data/mpra/k562_mpra_snps
python create_reporter_constructs.py media-3_oligos_snps_cleaned.tsv.gz ../data/mpra/k562_mpra_snps_2114 --procapnet

python create_reporter_constructs.py media-3_oligos_snps_cleaned.tsv.gz ../data/mpra/k562_mpra_snps --dinuc_shuffles 1
python create_reporter_constructs.py media-3_oligos_snps_cleaned.tsv.gz ../data/mpra/k562_mpra_snps_2114 --dinuc_shuffles 1 --procapnet
```

## Now predict on these sequences

We'll use the `predict_ensemble.py` script provided in the [original CLIPNET repo](https://github.com/Danko-Lab/clipnet/).

```bash
conda activate clipnet

cd /home2/ayh8/clipnet/ # Or wherever your CLIPNET path is.
# Switch paths to point to where your data and models are located.

# Example scripts for predicting with the model ensemble:
python predict_ensemble.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref.npz \
    --model_dir /home2/ayh8/clipnet_k562/models/clipnet_k562/ \
    --gpu 0
python predict_ensemble.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt.npz \
    --model_dir /home2/ayh8/clipnet_k562/models/clipnet_k562/ \
    --gpu 0

# To predict with the individual models:
for i in {1..9}; do
    python predict_individual_model.py \
        /home2/ayh8/clipnet_k562/models/clipnet_k562/fold_${i}.h5 \
        /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref.fa.gz \
        /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref_fold_${i}.h5 \
        --gpu 0;
done
for i in {1..9}; do
    python predict_individual_model.py \
        /home2/ayh8/clipnet_k562/models/clipnet_k562/fold_${i}.h5 \
        /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt.fa.gz \
        /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt_fold_${i}.h5 \
        --gpu 0;
done
```

## Predict using reference-trained CLIPNET

These models can be predicted with in the exact same fashion as above:

```bash
python predict_ensemble.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref_ref_model.h5 \
    --model_dir /home2/ayh8/clipnet_k562/models/clipnet_k562_reference/ \
    --gpu 0
python predict_ensemble.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt_ref_model.h5 \
    --model_dir /home2/ayh8/clipnet_k562/models/clipnet_k562_reference/ \
    --gpu 0

for i in {1..9}; do
    python predict_individual_model.py \
        /home2/ayh8/clipnet_k562/models/clipnet_k562_reference/fold_${i}.h5 \
        /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref.fa.gz \
        /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref_fold_${i}_ref_model.h5 \
        --gpu 0;
done
for i in {1..9}; do
    python predict_individual_model.py \
        /home2/ayh8/clipnet_k562/models/clipnet_k562_reference/fold_${i}.h5 \
        /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt.fa.gz \
        /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_alt_fold_${i}_ref_model.h5 \
        --gpu 0;
done
```

## Predict using ProCapNet

I wrote a short wrapper script to generate predictions w/ ProCapNet.

```bash
conda activate bpnet # environment into which you've installed the procapnet dependencies.

# cd to wherever you have this directory installed.
cd /home2/ayh8/clipnet_k562/mpra_benchmarks/
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
