"""
Calculates dataset parameters needed by clipnet for the K562 dataset.
"""

import json
import multiprocessing as mp
import os

import numpy as np

datadir = "/home2/ayh8/data/k562/k562_data_folds/"


def write_dataset_params(i):
    outdir = f"/home2/ayh8/k562_ensemble_models/f{i + 1}/"
    os.makedirs(outdir, exist_ok=True)

    test_folds = [i + 1]
    val_folds = [(i + 1) % 9 + 1]
    train_folds = [j for j in range(1, 10) if j not in test_folds + val_folds]
    print(train_folds, val_folds, test_folds)

    dataset_params = {
        "train_seq": [
            os.path.join(datadir, f"k562_sequence_{fold}.npz") for fold in train_folds
        ],
        "train_procap": [
            os.path.join(datadir, f"k562_procap_{fold}.npz") for fold in train_folds
        ],
        "val_seq": [
            os.path.join(datadir, f"k562_sequence_{fold}.npz") for fold in val_folds
        ],
        "val_procap": [
            os.path.join(datadir, f"k562_procap_{fold}.npz") for fold in val_folds
        ],
        "test_seq": [
            os.path.join(datadir, f"k562_sequence_{fold}.npz") for fold in test_folds
        ],
        "test_procap": [
            os.path.join(datadir, f"k562_procap_{fold}.npz") for fold in test_folds
        ],
    }

    dataset_params["n_train_folds"] = len(train_folds)
    dataset_params["n_val_folds"] = len(val_folds)
    dataset_params["n_test_folds"] = len(test_folds)

    # Calculate n_samples_per_chunk
    dataset_params["n_samples_per_train_fold"] = [
        np.load(f)["arr_0"].shape[0] for f in dataset_params["train_procap"]
    ]
    dataset_params["n_samples_per_val_fold"] = [
        np.load(f)["arr_0"].shape[0] for f in dataset_params["val_procap"]
    ]
    dataset_params["n_samples_per_test_fold"] = [
        np.load(f)["arr_0"].shape[0] for f in dataset_params["test_procap"]
    ]

    dataset_params["window_length"] = np.load(dataset_params["train_seq"][0])["arr_0"][
        0
    ].shape[0]

    dataset_params["pad"] = int(dataset_params["window_length"] / 4)
    dataset_params["output_length"] = int(
        2 * (dataset_params["window_length"] - 2 * dataset_params["pad"])
    )

    dataset_params["weight"] = 1 / 500

    output_fp = os.path.join(outdir, "dataset_params.json")

    with open(output_fp, "w") as handle:
        json.dump(dataset_params, handle, indent=4, sort_keys=True)


with mp.Pool(9) as p:
    p.map(write_dataset_params, range(9))
