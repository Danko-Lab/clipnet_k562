# Genomic benchmarks

Benchmark K562 models on predicting PRO-cap at genomic loci. To do this, we need to download K562 data

```bash
conda activate clipnet
# Requires that CLIPNET repo be cloned (https://github.com/Danko-Lab/clipnet) and the environment be installed.
cd ../../clipnet/

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
