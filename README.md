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
```

### Calculate CLIPNET_K562 attributions

```bash
conda activate clipnet
cd /home2/ayh8/clipnet/

python calculate_deepshap.py \
    ../clipnet_k562/data/k562_sequence_0.fa.gz \
    ../clipnet_k562/attributions/k562_deepshap_quantity_0.npz \
    ../clipnet_k562/attributions/k562_deepshap_quantity_0.npz
    --model_dir ../clipnet_k562/models/clipnet_k562 \
    --gpu 0 --mode 0
```
