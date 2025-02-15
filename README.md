# Analysis of transcription initiation models in K562

This repo contains code to analyse a couple different transcription initiation models in K562. Notably, we implement fine tuning to adapt CLIPNET models (trained in LCLs) to K562 and compare the ability of CLIPNET to predict MPRA SNP effects against ProCapNet, a transcription initiation model natively trained in K562, and Enformer, a multitask model trained on many K562 epigenetic tracks.

## Code dependencies

The python package requirements for CLIPNET K562 are listed in `requirements_tf.txt` and are identical to those used for the original CLIPNET project (https://github.com/Danko-Lab/clipnet/blob/main/requirements.txt).

ProCapNet (https://github.com/kundajelab/ProCapNet) and Enformer (https://github.com/lucidrains/enformer-pytorch) require Pytorch. The dependencies are listed in `requirements_pth.txt`. These should be installed in a separate environment from the TF environment used for CLIPNET. You will also need to install [PersonalBPNet](https://github.com/adamyhe/PersonalBPNet/), as it contains a number of utility functions and a copy of the ProCapNet class.

## Training data processing

Scripts to download/process K562 PRO-cap data for fine tuning are in `data_processing`. We note that that directory contains a number of pipelines for processing K562 data (including for PRO-seq models that have not yet been published). For more details on the minimal necessary scripts to reproduce the K562 initiation models, please consult the README in that directory.

## Fine tuning scripts

Scripts to fine tune CLIPNET models to K562 are located in `clipnet_ft`. These scripts will presume that the training data have been downloaded from Zenodo or preprocessed per protocols in `data_processing/`

## Download pretrained models

### CLIPNET K562

These models are compatible with the scripts in the original CLIPNET repo (simply use `-m models/clipnet_k562/` to specify the directory into which these models have been downloaded). For more details, see `clipnet_ft`.

```bash
mkdir -p models/clipnet_k562/
for fold in {1..9}; do
    wget https://zenodo.org/records/11196189/files/fold_${fold}.h5 -P models/clipnet_k562/;
done
```

### ProCapNet K562

NOTE: ProCapNet and CLIPNET use different libraries, so you should install the environment for each model separately.

```bash
mkdir -p models/procapnet_k562/
wget https://www.encodeproject.org/files/ENCFF976FHE/@@download/ENCFF976FHE.tar.gz
tar -xvf ENCFF976FHE.tar.gz -C models/procapnet_k562/
rm ENCFF976FHE.tar.gz
```

## Benchmarking

Scripts to benchmark models on PRO-cap prediction at genomic loci are in `genomic_benchmarks`.

Scripts to benchmark models on MPRA variant effect prediction are in `mpra_benchmarks`.
