# Analysis of transcription initiation models in K562

Transfer learning of CLIPNET (trained in LCLs) to K562.

Scripts to download/process K562 PRO-cap data for transfer learning are in `data_processing`. Scripts to transfer learn CLIPNET models to K562 are in `clipnet_transfer_learning`. Scripts to benchmark models across loci are in `clipnet_transfer_learning`. Scripts to benchmark models on MPRA variant prediction are in `siraj_mpra`.

## Download pretrained models

### CLIPNET K562

```bash
for fold in {1..9};
do wget https://zenodo.org/records/11196189/files/fold_${fold}.h5 -P models/clipnet_k562/;
done
```

### ProCapNet K562

```bash
# NOTE: ProCapNet and CLIPNET use different libraries, so you should install the environment for each model separately.
wget https://www.encodeproject.org/files/ENCFF976FHE/@@download/ENCFF976FHE.tar.gz
tar -xvf ENCFF976FHE.tar.gz -C models/procapnet_k562/
rm ENCFF976FHE.tar.gz
```
