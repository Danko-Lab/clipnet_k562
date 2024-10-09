# QTL analysis of MPRA data from Siraj et al. 2024

This repository contains the code used to perform the QTL analysis of MPRA data from [Siraj et al. 2024](http://dx.doi.org/10.1101/2024.05.05.592437). The data used in this analysis was downloaded from the supplementary material of the paper.
$$
Relevant addgene sequences:

https://www.addgene.org/109035/sequences/

## Process SI data

Use `process_SI.ipynb`.

## Generate reporter constructs

Run `create_reporter_constructs.py`. This script will generate the reporter constructs used in the MPRA experiment (may be worth confirming with the original authors), and random sequences to be use to potentially marginalize out minP activity.

```bash
python create_reporter_constructs_for_fine_tuning.py ../siraj_mpra/media-3_oligos_snps.tsv.gz ../../data/k562_mpra/k562_mpra_snps_ft
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
python predict_ensemble.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_shuffle1.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_shuffle1.h5 \
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
```
