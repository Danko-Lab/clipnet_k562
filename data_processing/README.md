# Fixed test window pipelines

## Summary

The goal of this pipeline is to generate a fixed set of windows benchmark the performance of each of the CLIPNET replicates. Peaks were selected such that they are present in at least 60 libraries. These pipelines should be run after the main `lcl_snakemake_pipelines` have been run to generate consensus sequences, peak calls, etc.

## Get windows

```bash
cd ./get_window
bash make.sh
cd ../
```

## Process PRO-cap

```bash
cd ./procap
bash make.sh
cd ../
```

## Process sequence

```bash
cd ./sequence
bash make.sh
cd ../
```

## Predict and run performance metrics

```bash
cd /home2/ayh8/clipnet/
conda activate clipnet

python predict_ensemble.py \
    /home2/ayh8/data/k562/k562_data_folds/k562_sequence_0.fna.gz \
    /home2/ayh8/predictions/k562/k562_prediction_0.h5 \
    --model_dir /home2/ayh8/k562_ensemble_models/ \
    --gpu
python calculate_performance_metrics.py \
    /home2/ayh8/predictions/k562/k562_prediction_0.h5 \
    /home2/ayh8/data/k562/k562_data_folds/k562_procap_0.npz \
    /home2/ayh8/predictions/k562/k562_0erformance_0.h5
```
