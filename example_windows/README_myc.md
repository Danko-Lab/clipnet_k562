# Prediction at myc promoter

## Extract window

```bash
# same center as Cochran et al.
twoBitToFa /fs/cbsubscb17/storage/data/hg38/hg38.2bit:chr8:127735632-127736632 /home2/ayh8/data/myc_promoter/myc_promoter_0.fa

# slightly more downstream
twoBitToFa /fs/cbsubscb17/storage/data/hg38/hg38.2bit:chr8:127735655-127736655 /home2/ayh8/data/myc_promoter/myc_promoter_1.fa
```

## Predict

```bash
cd /home2/ayh8/clipnet
conda activate clipnet

i=9
# LCL model
python predict_ensemble.py \
    /home2/ayh8/data/myc_promoter/myc_promoter_0.fa \
    /home2/ayh8/data/myc_promoter/myc_promoter_0_ensemble.h5
python predict_ensemble.py \
    /home2/ayh8/data/myc_promoter/myc_promoter_1.fa \
    /home2/ayh8/data/myc_promoter/myc_promoter_1_ensemble.h5
python predict_individual_model.py \
    ensemble_models/fold_${i}.h5 \
    /home2/ayh8/data/myc_promoter/myc_promoter_0.fa \
    /home2/ayh8/data/myc_promoter/myc_promoter_0_fold_${i}.h5
python predict_individual_model.py \
    ensemble_models/fold_${i}.h5 \
    /home2/ayh8/data/myc_promoter/myc_promoter_1.fa \
    /home2/ayh8/data/myc_promoter/myc_promoter_1_fold_${i}.h5

# CD4 model
python predict_ensemble.py \
    /home2/ayh8/data/myc_promoter/myc_promoter_0.fa \
    /home2/ayh8/data/myc_promoter/myc_promoter_0_ensemble_cd4.h5 \
    --model_dir ../cd4_ensemble_models_0406/
python predict_ensemble.py \
    /home2/ayh8/data/myc_promoter/myc_promoter_1.fa \
    /home2/ayh8/data/myc_promoter/myc_promoter_1_ensemble_cd4.h5 \
    --model_dir ../cd4_ensemble_models_0406/
python predict_individual_model.py \
    ../cd4_ensemble_models_0406/fold_${i}.h5 \
    /home2/ayh8/data/myc_promoter/myc_promoter_0.fa \
    /home2/ayh8/data/myc_promoter/myc_promoter_0_fold_${i}_cd4.h5
python predict_individual_model.py \
    ../cd4_ensemble_models_0406/fold_${i}.h5 \
    /home2/ayh8/data/myc_promoter/myc_promoter_1.fa \
    /home2/ayh8/data/myc_promoter/myc_promoter_1_fold_${i}_cd4.h5

# K562 model
python predict_ensemble.py \
    /home2/ayh8/data/myc_promoter/myc_promoter_0.fa \
    /home2/ayh8/data/myc_promoter/myc_promoter_0_ensemble_k562.h5 \
    --model_dir ../clipnet_k562/models/clipnet_k562/
python predict_ensemble.py \
    /home2/ayh8/data/myc_promoter/myc_promoter_1.fa \
    /home2/ayh8/data/myc_promoter/myc_promoter_1_ensemble_k562.h5 \
    --model_dir ../clipnet_k562/models/clipnet_k562/
python predict_individual_model.py \
    ../clipnet_k562/models/clipnet_k562/fold_${i}.h5 \
    /home2/ayh8/data/myc_promoter/myc_promoter_0.fa \
    /home2/ayh8/data/myc_promoter/myc_promoter_0_fold_${i}_k562.h5
python predict_individual_model.py \
    ../clipnet_k562/models/clipnet_k562/fold_${i}.h5 \
    /home2/ayh8/data/myc_promoter/myc_promoter_1.fa \
    /home2/ayh8/data/myc_promoter/myc_promoter_1_fold_${i}_k562.h5
```
