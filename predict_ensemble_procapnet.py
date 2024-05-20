import procapnet
import numpy as np
import h5py

def predict_ensemble_procapnet(X, model_dir, batch_size=64, logits=False):