# CLIPNET K562

Transfer learning of CLIPNET (trained in LCLs) to K562

## Download pretrained models

### CLIPNET_K562

```bash
```

### ProCapNet_K562

```bash
```

### Calculate CLIPNET_K562 predictions

```bash
conda activate clipnet
cd /home2/ayh8/clipnet/

mkdir -p ../clipnet_k562/predictions/
python predict_ensemble.py \
    ../clipnet_k562/data/k562_data_folds/k562_sequence_0.fa.gz \
    ../clipnet_k562/predictions/k562_predictions_0.h5 \
    --model_dir ../clipnet_k562/models/clipnet_k562 \
    --gpu 0

mkdir -p ../clipnet_k562/predictions/ensemble_test/
for i in {1..9}; do
    python predict_individual_model.py \
        ../clipnet_k562/models/clipnet_k562/fold_${i}.h5 \
        ../clipnet_k562/data/k562_data_folds/k562_sequence_0.fa.gz \
        ../clipnet_k562/predictions/ensemble_test/k562_predictions_${i}.h5 \
        --gpu 0;
done

mkdir -p ../clipnet_k562/predictions/individual_test/
for i in {1..9}; do
    python predict_individual_model.py \
        ../clipnet_k562/models/clipnet_k562/fold_${i}.h5 \
        ../clipnet_k562/data/k562_data_folds/k562_sequence_${i}.fa.gz \
        ../clipnet_k562/predictions/individual_test/k562_predictions_${i}.h5 \
        --gpu 0;
done

python calculate_performance_metrics.py \
    ../clipnet_k562/predictions/k562_predictions_0.h5 \
    ../clipnet_k562/data/k562_data_folds/k562_procap_0.npz \
    ../clipnet_k562/predictions/k562_performance_0.h5

for i in {1..9}; do
    python calculate_performance_metrics.py \
        ../clipnet_k562/predictions/ensemble_test/k562_predictions_${i}.h5 \
        ../clipnet_k562/data/k562_data_folds/k562_procap_0.npz \
        ../clipnet_k562/predictions/ensemble_test/k562_performance_${i}.h5
done

for i in {1..9}; do
    python calculate_performance_metrics.py \
        ../clipnet_k562/predictions/individual_test/k562_predictions_${i}.h5 \
        ../clipnet_k562/data/k562_data_folds/k562_procap_${i}.npz \
        ../clipnet_k562/predictions/individual_test/k562_performance_${i}.h5
done
```

### Calculate CLIPNET_K562 attributions

```bash
conda activate clipnet
cd /home2/ayh8/clipnet/

mode=profile
python calculate_deepshap.py \
    ../clipnet_k562/data/k562_data_folds/k562_centered.fa.gz \
    ../clipnet_k562/attributions/k562_deepshap_${mode}.npz \
    ../clipnet_k562/attributions/k562_seqs_onehot.npz \
    --model_dir ../clipnet_k562/models/clipnet_k562 \
    --gpu 0 --mode $mode
```

### Calculate CLIPNET_K562 TF-MoDISco

```bash
conda activate modisco
mode=profile

time modisco motifs -s attributions/k562_seqs_onehot.npz -a attributions/k562_deepshap_${mode}.npz -n 1000000 -l 50 -v -o attributions/k562_deepshap_${mode}_modisco.h5
time modisco report -i attributions/k562_deepshap_${mode}_modisco.h5 -o attributions/k562_deepshap_${mode}_modisco/ -m /home2/ayh8/data/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
```
