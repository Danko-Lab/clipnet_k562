# Transfer learning CLIPNET to K562

This directory contains scripts for transfer learning the CLIPNET model (trained in LCLs) to K562 cells.

First, download and process data as described in `data_processing/README.md`.

Then the following scripts can be used to train the model:

`calculate_fold_params.py` calculates the fold parameters for the model (how many examples/epoch, sets up all the file names correctly, etc).

`transfer_learn_k562.py` trains the model.

`transfer_learn_k562_reference.py` transfer learns the reference-trained CLIPNET model.

`utils.py` contains utility functions for training the model.

Other scripts are for training a PRO-seq model in K562 cells (still in-progress).

Unless you want to replicate the results of the study, you should not need to run these scripts. The pre-trained model weights are available on Zenodo at https://zenodo.org/records/11196189.