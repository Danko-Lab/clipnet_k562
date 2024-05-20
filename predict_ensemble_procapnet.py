import procapnet
import numpy as np
import h5py
import glob
import utils

paths = glob.glob("models/procapnet_k562/*.model")
models = [procapnet.Model(p) for p in paths]

ref_seqs = utils.get_twohot_fasta_sequences("data/mpra/k562_mpra_snps_ref.fa.gz") / 2
alt_seqs = utils.get_twohot_fasta_sequences("data/mpra/k562_mpra_snps_alt.fa.gz") / 2

