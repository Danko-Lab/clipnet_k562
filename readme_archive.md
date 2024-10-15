### Calculate CLIPNET_K562 predictions

```bash
conda activate clipnet
# Requires that CLIPNET repo be cloned (https://github.com/Danko-Lab/clipnet)
cd ../clipnet/

mkdir -p ../predictions/k562/
python predict_ensemble.py \
    ../data/k562_data_folds/k562_sequence_0.fna.gz \
    ../predictions/k562/k562_predictions_0.h5 \
    --model_dir ../clipnet_k562/models/clipnet_k562 \
    --gpu 0

mkdir -p ../predictions/k562/individual_test/
for i in {1..9}; do
    python predict_individual_model.py \
        ../clipnet_k562/models/clipnet_k562/fold_${i}.h5 \
        ../data/k562_data_folds/k562_sequence_${i}.fna.gz \
        ../predictions/k562/individual_test/k562_predictions_${i}.h5 \
        --gpu 0;
done

python calculate_performance_metrics.py \
    ../predictions/k562/k562_predictions_0.h5 \
    ../data/k562_data_folds/k562_procap_0.npz \
    ../predictions/k562/k562_performance_0.h5

for i in {1..9}; do
    python calculate_performance_metrics.py \
        ../predictions/k562/individual_test/k562_predictions_${i}.h5 \
        ../data/k562_data_folds/k562_procap_${i}.npz \
        ../predictions/k562/individual_test/k562_performance_${i}.h5
done
```

### Calculate CLIPNET_K562 PRO-seq predictions

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

### Calculate CLIPNET_K562 attributions

```bash
conda activate clipnet
cd /home2/ayh8/clipnet/

mode=profile
python calculate_deepshap.py \
    ../data/k562_data_folds/k562_data_folds/k562_centered.fna.gz \
    ../clipnet_k562/attributions/k562_deepshap_${mode}.npz \
    ../clipnet_k562/attributions/k562_seqs_onehot.npz \
    --model_dir ../clipnet_k562/models/clipnet_k562 \
    --gpu 0 --mode $mode
```

### Calculate CLIPNET_K562 PRO-seq attributions

```bash
conda activate clipnet
cd /home2/ayh8/clipnet/

mode=profile
python calculate_deepshap.py \
    ../data/k562_data_folds/k562_data_folds/k562_centered.fna.gz \
    ../clipnet_k562/attributions/k562_deepshap_${mode}.npz \
    ../clipnet_k562/attributions/k562_seqs_onehot.npz \
    --model_dir ../clipnet_k562/models/clipnet_k562 \
    --gpu 0 --mode $mode
```

### Calculate CLIPNET_K562 TF-MoDISco

```bash
conda activate modisco
mode=profile

time modisco motifs \
    -s attributions/k562_seqs_onehot.npz \
    -a attributions/k562_deepshap_${mode}.npz \
    -o attributions/k562_deepshap_${mode}_modisco.h5 \
    -n 1000000 -l 50 -v 
time modisco report \
    -i attributions/k562_deepshap_${mode}_modisco.h5 \
    -o attributions/k562_deepshap_${mode}_modisco/ \
    -m /home2/ayh8/data/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
```
