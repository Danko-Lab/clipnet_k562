# Fine tuning CLIPNET to K562

This directory contains scripts for fine tuning the CLIPNET model (trained in LCLs) to K562 data.

Unless you want to replicate the results of the study, you should not need to run these scripts. The pre-trained model weights are available on Zenodo at https://zenodo.org/records/11196189.

## Fine tuning

First, download and process data as described in `data_processing/README.md`.

Then the following scripts can be used to train the model:

`calculate_fold_params.py` calculates the fold parameters for the model (how many examples/epoch, sets up all the file names correctly, etc).

`ft_k562.py` fine tunes the model to K562.

`ft_k562_reference.py` fine tunes the reference-trained CLIPNET model.

To do the fine tuning, download and unpack pre-trained LCL models. Then run `calculate_fold_params.py` to calculate the fold parameters. Finally, run `ft_k562.py` to train the model.

```bash
# model download
mkdir -p ../models/clipnet_lcl/
for fold in {1..9}; do
    wget https://zenodo.org/records/10408623/files/fold_${fold}.h5 -P ../models/clipnet_lcl/;
done

# DATADIR = where processed K562 data is stored.
mkdir -p ../models/clipnet_k562/
python calculate_fold_params.py $DATADIR ../models/clipnet_k562/
GPU=0
for i in {1..9}; do python ft_k562.py $i $GPU; done
```

We can also fine tune from the reference-trained LCL CLIPNET models:

```bash
# model download
mkdir -p ../models/clipnet_lcl_reference/
wget https://zenodo.org/records/14037356/files/reference_models.tar
tar -xvf reference_models.tar -C ../models/clipnet_lcl_reference/ --strip-components=1
rm reference_models.tar

# DATADIR = where processed K562 data is stored.
mkdir -p ../models/clipnet_k562/
python calculate_fold_params.py $DATADIR ../models/clipnet_k562_reference/
GPU=0
for i in {1..9}; do python ft_k562_reference.py $i $GPU; done
```

## K562 model downloads

The fine tuned models that we use in our paper are deposited at Zenodo. These can be downloaded as follows (obviously if you run this you will not need to do your own fine tuning):

For the K562 models fine tuned from the personal CLIPNET LCL models (recommended for general use):

```bash
mkdir -p ../models/clipnet_k562/
for fold in {1..9}; do
    wget https://zenodo.org/records/11196189/files/fold_${fold}.h5 -P ../models/clipnet_k562/;
done
```

For the K562 models fine tuned from the reference-trained CLIPNET LCL models (recommended only if you really want no genetic variation data leakage on genome-wide VEP):

```bash
mkdir -p ../models/clipnet_k562_reference/
wget https://zenodo.org/records/14037356/files/clipnet_k562_reference.tar
tar -xvf clipnet_k562_reference.tar -C ../models/clipnet_k562_test --strip-components=1
rm clipnet_k562_reference.tar
```
