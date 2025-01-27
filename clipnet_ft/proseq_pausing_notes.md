# PRO-seq and pausing index model notes

These describe unpublished work that is not well documented (sorry), but I'm keeping this here for convenience.

We fine tune the PRO-cap model to PRO-seq data:

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
