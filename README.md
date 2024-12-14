# Analysis of transcription initiation models in K562

This repo contains code to analyse a couple different transcription initiation models in K562. Notably, we implement transfer learning to adapt CLIPNET models CLIPNET (trained in LCLs) to K562 and compare the ability of CLIPNET to predict MPRA SNP effects against ProCapNet, a transcription initiation model natively trained in K562, and Enformer, a multitask model trained on many K562 epigenetic tracks.

## Code dependencies

The python package requirements for CLIPNET K562 are listed in `requirements_tf.txt` and are identical to those used for the original CLIPNET project (https://github.com/Danko-Lab/clipnet/blob/main/requirements.txt).

ProCapNet (https://github.com/kundajelab/ProCapNet) and Enformer (https://github.com/lucidrains/enformer-pytorch) require Pytorch. The dependencies are listed in `requirements_pth.txt`. These should be installed in a separate environment from the TF environment used for CLIPNET.

## Training data processing

Scripts to download/process K562 PRO-cap data for transfer learning are in `data_processing`. We note that that directory contains a number of pipelines for processing K562 data (including for PRO-seq models that have not yet been published). For more details on the minimal necessary scripts to reproduce the K562 initiation models, please consult the README in that directory.

## Transfer learning scripts

Scripts to transfer learn CLIPNET models to K562 and benchmarking their performance at predicting initiation across genomic loci are located in `clipnet_transfer_learning`.

## Download pretrained models

### CLIPNET K562

```bash
for fold in {1..9};
do wget https://zenodo.org/records/11196189/files/fold_${fold}.h5 -P models/clipnet_k562/;
done
```

### ProCapNet K562

```bash
# NOTE: ProCapNet and CLIPNET use different libraries,
# so you should install the environment for each model separately.
wget https://www.encodeproject.org/files/ENCFF976FHE/@@download/ENCFF976FHE.tar.gz
tar -xvf ENCFF976FHE.tar.gz -C models/procapnet_k562/
rm ENCFF976FHE.tar.gz
```

Scripts to benchmark models on MPRA variant prediction are in `siraj_mpra`.