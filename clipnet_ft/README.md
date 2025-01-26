# Fine tuning CLIPNET to K562

This directory contains scripts for fine tuning the CLIPNET model (trained in LCLs) to K562 data.

Unless you want to replicate the results of the study, you should not need to run these scripts. The pre-trained model weights are available on Zenodo at https://zenodo.org/records/11196189.

## Fine tuning

First, download and process data as described in `data_processing/README.md`.

Then the following scripts can be used to train the model:

`calculate_fold_params.py` calculates the fold parameters for the model (how many examples/epoch, sets up all the file names correctly, etc).

`ft_k562.py` fine tunes the model to K562.

`ft_k562_reference.py` fine tunes the reference-trained CLIPNET model.

To do the fine tuning, first run `calculate_fold_params.py` to calculate the fold parameters. Then run `ft_k562.py` to train the model:

```bash
#DATADIR = where data is stored.
#OUTDIR = where models will be saved
python calculate_fold_params.py $DATADIR $OUTDIR
GPU=0
for i in {1..9}; do python ft_k562.py $i $GPU; done
```

We then fine tune the PRO-cap model to PRO-seq data (unpublished work, not extensively documented here):

```bash
python calculate_fold_params_proseq.py $DATADIR $OUTDIR
GPU=0
for i in {1..9}; do python ft_k562_proseq.py $i $GPU; done
```

And then the PRO-seq models to predict pausing index at promoters (also not documented here). First, calculate pausing indices and get promoter sequences using pipelines in `data_processing/pausing_index` and `data_processing/pausing_index_sequence`.

```bash
for i in {1..9}; do python ft_k562_pausing.py $i; done
```

Calculate attributions for the pausing model:

```bash
for i in {1..9}; do python pausing_deepshap.py ../../data/pausing_index/k562_pausing_index_centered.fa.gz ../../data/pausing_index/k562_pausing_index_centered_deepshap_${i}.npz --ohe_seq_fp ../../data/k562_pausing_index_centered_ohe.npz --model_fp ../models/clipnet_k562_pausing/fold_${i}.h5; done
```

Average DeepSHAP files, then calculate MODISCO:

```bash
cd ../../data/pausing_index/
time modisco motifs -s k562_pausing_index_centered_ohe.npz -a k562_pausing_index_centered_deepshap.npz -n 1000000 -l 50 -v -o k562_pausing_index_centered_modisco.h5
time modisco report -i k562_pausing_index_centered_modisco.h5 -o k562_pausing_index_centered_modisco/ -m /home2/ayh8/data/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
```
