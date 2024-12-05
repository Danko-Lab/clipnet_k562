# Transfer learning CLIPNET to K562

This directory contains scripts for transfer learning the CLIPNET model (trained in LCLs) to K562 data.

Unless you want to replicate the results of the study, you should not need to run these scripts. The pre-trained model weights are available on Zenodo at https://zenodo.org/records/11196189.

## Transfer learning

First, download and process data as described in `data_processing/README.md`.

Then the following scripts can be used to train the model:

`calculate_fold_params.py` calculates the fold parameters for the model (how many examples/epoch, sets up all the file names correctly, etc).

`transfer_learn_k562.py` trains the model.

`transfer_learn_k562_reference.py` transfer learns the reference-trained CLIPNET model.

To do the transfer learning, first run `calculate_fold_params.py` to calculate the fold parameters. Then run `transfer_learn_k562.py` to train the model:

```bash
python calculate_fold_params.py $DATADIR $OUTDIR
GPU=0
for i in {1..9}; do python transfer_learn_k562.py $i $GPU; done
```
