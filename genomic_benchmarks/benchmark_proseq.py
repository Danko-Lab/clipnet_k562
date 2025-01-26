import h5py
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.stats import pearsonr, spearmanr


def average_by_modulo(arr, m, ax=0):
    # Get the shape of the array (mx3, n)
    total_rows, n_cols = arr.shape
    # Calculate the number of groups (residue classes 0, 1, ..., m-1)
    groups = m
    # Initialize an array to hold the averages for each residue class
    result = np.zeros((m, n_cols))
    # Compute averages by residue class modulo m
    for i in range(m):
        result[i] = np.mean(arr[i::m], axis=ax)  # Select rows where index % m == i
    return result


prediction = h5py.File("../../predictions/k562/k562_proseq_nonintragenic.h5")
pred = (
    prediction["track"][:]
    / prediction["track"][:].sum(axis=1, keepdims=True)
    * prediction["quantity"][:]
)
observed = average_by_modulo(
    np.load("../../data/k562/k562_data_folds/k562_proseq_0.npz")["arr_0"],
    prediction["track"].shape[0],
)

track_corr = pd.DataFrame(observed[:, np.r_[250:750, 1250:1750]]).corrwith(
    pd.DataFrame(pred), axis=1
)
jsd = jensenshannon(observed[:, np.r_[250:750, 1250:1750]], pred)
pearsonr(
    np.log(observed[:, np.r_[250:750, 1250:1750]].sum(axis=1) + 1),
    np.log(prediction["quantity"][:, 0] + 1),
)
spearmanr(
    observed[:, np.r_[250:750, 1250:1750]].sum(axis=1), prediction["quantity"][:, 0]
)
