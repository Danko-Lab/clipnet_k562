"""
Calculates dataset parameters needed by clipnet. Supply a path to the processed data
and an output directory. This script will write json files to the output directory
for use in model training.
"""

import argparse
import json
import os

import numpy as np


def write_dataset_params(i, datadir, outdir):
    outdir = f"{outdir}/f{i + 1}/"
    os.makedirs(outdir, exist_ok=True)

    test_folds = [i + 1]
    val_folds = [(i + 1) % 9 + 1]
    train_folds = [j for j in range(1, 10) if j not in test_folds + val_folds]
    print(train_folds, val_folds, test_folds)

    dataset_params = {
        "train_seq": [
            os.path.join(datadir, f"k562_sequence_nonintragenic_{fold}.npz")
            for fold in train_folds
        ],
        "train_proseq": [
            os.path.join(datadir, f"k562_proseq_{fold}.npz") for fold in train_folds
        ],
        "val_seq": [
            os.path.join(datadir, f"k562_sequence_nonintragenic_{fold}.npz")
            for fold in val_folds
        ],
        "val_proseq": [
            os.path.join(datadir, f"k562_proseq_{fold}.npz") for fold in val_folds
        ],
        "test_seq": [
            os.path.join(datadir, f"k562_sequence_nonintragenic_{fold}.npz")
            for fold in test_folds
        ],
        "test_proseq": [
            os.path.join(datadir, f"k562_proseq_{fold}.npz") for fold in test_folds
        ],
    }

    dataset_params["n_train_folds"] = len(train_folds)
    dataset_params["n_val_folds"] = len(val_folds)
    dataset_params["n_test_folds"] = len(test_folds)

    # Calculate n_samples_per_chunk
    dataset_params["n_samples_per_train_fold"] = [
        np.load(f)["arr_0"].shape[0] for f in dataset_params["train_proseq"]
    ]
    dataset_params["n_samples_per_val_fold"] = [
        np.load(f)["arr_0"].shape[0] for f in dataset_params["val_proseq"]
    ]
    dataset_params["n_samples_per_test_fold"] = [
        np.load(f)["arr_0"].shape[0] for f in dataset_params["test_proseq"]
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("datadir", type=str, help="directory containing data")
    parser.add_argument(
        "outdir",
        type=str,
        help="directory to save dataset params to (where models will be saved)",
    )
    parser.add_argument(
        "--threads", type=int, default=9, help="number of threads to use"
    )
    args = parser.parse_args()
    if args.threads == 1:
        for i in range(9):
            write_dataset_params(i, args.datadir, args.outdir)
    elif args.threads > 1:
        import itertools
        import multiprocessing as mp

        with mp.Pool(9) as p:
            p.starmap(
                write_dataset_params,
                zip(
                    range(9),
                    itertools.repeat(args.datadir),
                    itertools.repeat(args.outdir),
                ),
            )
    else:
        raise ValueError("threads must be >= 1")


if __name__ == "__main__":
    main()
