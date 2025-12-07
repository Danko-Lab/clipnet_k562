# ARCHIVED SCRIPTS

This md contains a bunch of scripts for generating data that didn't make into a paper.

## Calculate CLIPNET_K562 PRO-seq predictions

```bash
conda activate clipnet
cd /home2/ayh8/clipnet/

mkdir -p ../predictions/k562/
python predict_ensemble.py \
    ../data/k562_data_folds/k562_sequence_nonintragenic_0.fna.gz \
    ../predictions/k562/k562_proseq_predictions_0.h5 \
    --model_dir ../clipnet_k562/models/clipnet_k562_proseq \
    --gpu 0

python calculate_performance_metrics.py \
    ../predictions/k562/k562_proseq_predictions_0.h5 \
    ../data/k562_data_folds/k562_proseq_0.npz \
    ../predictions/k562/k562_proseq_performance_metrics_0.h5
```

## Calculate CLIPNET_K562 attributions

```bash
conda activate clipnet
cd /home2/ayh8/clipnet_k562/

attr=profile
clipnet attribute \
    -f ../data/k562/k562_centered.fa.gz \
    -o ../attr/k562/k562_procap_deepshap_${attr}.npz \
    -m models/clipnet_k562 \
    -a $attr -c -v
```

## Calculate CLIPNET_K562 PRO-seq attributions

```bash
conda activate clipnet
cd /home2/ayh8/clipnet_k562/

attr=profile
clipnet attribute \
    -f ../data/k562/k562_centered.fa.gz \
    -o ../attr/k562/k562_proseq_deepshap_${attr}.npz \
    -m models/clipnet_k562_proseq \
    -a $attr -c -v
```

## Calculate CLIPNET_K562 TF-MoDISco

```bash
conda activate bpnet
cd /home2/ayh8/attr/k562/

attr=profile
time modisco motifs \
    -s attributions/k562_seqs_onehot.npz \
    -a attributions/k562_proseq_deepshap_${attr}.npz \
    -o attributions/k562_proseq_deepshap_${attr}_modisco.h5 \
    -n 1000000 -l 50 -v 
time modisco report \
    -i attributions/k562_proseq_deepshap_${attr}_modisco.h5 \
    -o attributions/k562_proseq_deepshap_${attr}_modisco/ \
    -m ../data/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt
```
