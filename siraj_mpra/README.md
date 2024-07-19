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

## Predict using ProCapNet

```bash
conda activate procapnet

python predict_ensemble_procapnet.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_2114_ref.fa.gz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_2114_ref.h5 \
    --model_dir /home2/ayh8/clipnet_k562/models/procapnet_k562/ \
    --gpu 0
```

## Calculate performance metrics

```bash
python calculate_performance_metrics.py \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref.h5 \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref_procapnet.npz \
    /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_ref_performance.h5
```

## Split out data fold 0

```python
import pandas as pd
import pyfastx
import tqdm

data_folds = pd.read_csv("clipnet_data_fold_assignments.csv")
df = pd.read_csv("../data/mpra/k562_allelic_mpra_snps.csv.gz")
ref_fasta = pyfastx.Fasta("../data/mpra/k562_mpra_snps_ref.fa.gz")
alt_fasta = pyfastx.Fasta("../data/mpra/k562_mpra_snps_alt.fa.gz")

seq_df = pd.DataFrame(
    {
        "variant": [x.name.split("_")[0] for x in ref_fasta],
        "ref_seq": [x.seq for x in ref_fasta],
        "alt_seq": [x.seq for x in alt_fasta],
    }
)

test_df = df[df["fold"] == 0]
merged_df = seq_df.merge(test_df, on="variant")

test_ref_fasta_fp = "../data/mpra/k562_mpra_snps_ref_fold0.fa"
with open(test_ref_fasta_fp, "w") as f:
    for i, row in tqdm.tqdm(merged_df.iterrows(), total=len(merged_df)):
        f.write(f">{row['variant']}_ref\n{row['ref_seq']}\n")


test_alt_fasta_fp = "../data/mpra/k562_mpra_snps_alt_fold0.fa"
with open(test_alt_fasta_fp, "w") as f:
    for i, row in tqdm.tqdm(merged_df.iterrows(), total=len(merged_df)):
        f.write(f">{row['variant']}_alt\n{row['alt_seq']}\n")

test_df.to_csv("../data/mpra/k562_mpra_snps_fold0.csv.gz", index=False)
```

## Run DeepSHAP

```bash
conda activate clipnet

cd /home2/ayh8/clipnet/
allele=alt
mode=profile
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
allele=alt
mode=quantity
for fold in {0..6}; do 
    python calculate_deepshap_procapnet.py \
        /home2/ayh8/clipnet_k562/models/procapnet_k562/fold_${fold}/ \
        /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_${allele}_fold0.fa \
        /home2/ayh8/clipnet_k562/data/mpra/k562_mpra_snps_${allele}_fold0_procapnet_deepshap_${mode}.npz \
        --mode $mode \
        --gpu 1;
done
```
