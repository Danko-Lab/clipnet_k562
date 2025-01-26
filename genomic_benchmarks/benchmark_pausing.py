import logging
import os

import numpy as np
import pandas as pd
import tqdm
import utils
from scipy.stats import pearsonr, spearmanr

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "4"
logging.getLogger("tensorflow").setLevel(logging.FATAL)
import tensorflow as tf

# Set VRAM usage to growth:
gpus = tf.config.experimental.list_physical_devices("GPU")
for gpu in gpus:
    tf.config.experimental.set_memory_growth(gpu, True)

# Load pausing indices
pausing_index = pd.read_csv(
    "../../data/pausing_index/k562_pausing_index_G156.bed",
    header=None,
    sep="\t",
    names=("chrom", "start", "stop", "name", "y", "strand"),
)

# Load sequences
sequences = utils.get_twohot_fasta_sequences(
    "../../data/pausing_index/k562_pausing_index_centered.fa.gz"
)

# Extract test set
test_folds = [0]
chrom_splits = pd.read_csv("../data_processing/data_fold_assignments.csv")
test_chroms = [
    row["chrom"] for i, row in chrom_splits.iterrows() if row["fold"] in test_folds
]

test_idx = pausing_index.chrom.isin(test_chroms)
X = sequences[test_idx]
y = pausing_index[test_idx].y.to_numpy()

# Load models
models = [
    tf.keras.models.load_model(
        f"../models/clipnet_k562_pausing/fold_{i}.h5", compile=False
    )
    for i in tqdm.trange(1, 10, desc="Loading models")
]


# Predict
def predict(model, X, batch_size=64, silence=False):
    y_predict_handle = [
        model.predict(X[i : i + batch_size, :, :], verbose=0)
        for i in tqdm.tqdm(
            range(0, X.shape[0], batch_size),
            desc=f"Predicting in batches of {batch_size}",
            disable=silence,
        )
    ]
    y_predict = np.concatenate([chunk for chunk in y_predict_handle])
    return y_predict


y_pred = [predict(model, X) for model in models]

for p in y_pred:
    print(
        pearsonr(np.log(p.squeeze() + 1e-6), np.log(y + 1e-6)),
        spearmanr(p.squeeze(), y),
    )

y_pred_mean = np.stack(y_pred, axis=0).mean(axis=0)
print(
    pearsonr(np.log(y_pred_mean.squeeze() + 1e-6), np.log(y + 1e-6)),
    spearmanr(p.squeeze(), y),
)
np.savez_compressed(
    "../../data/pausing_index/k562_pausing_index_centered_prediction.npz",
    y_pred=y_pred_mean,
    y=y,
)

# PearsonRResult(statistic=0.5766288425062598, pvalue=0.0)
# SignificanceResult(statistic=0.6673798094684796, pvalue=0.0)
