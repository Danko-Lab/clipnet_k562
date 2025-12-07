"""
This script benchmarks the ProCapNet model on the K562 ProCap dataset.
"""

import argparse
import glob

import numpy as np
import pyfastx
import torch
import tqdm
from bpnetlite.attribute import deep_lift_shap
from bpnetlite.bpnet import CountWrapper, ProfileWrapper
from personal_bpnet import procapnet

# from tangermeme.deep_lift_shap import deep_lift_shap
from tangermeme.utils import one_hot_encode


def save_deepshap_results(onehot_seqs, scores, scores_path, onehot_seqs_path=None):
    if len(onehot_seqs.shape) != 3 or onehot_seqs.shape[1] != 4:
        raise ValueError(
            f"{onehot_seqs.shape} is incorrect shape for one-hot encoded sequences. "
            + "Expected (n_seqs, 4, seq_len)."
        )
    if len(scores.shape) != 3 or scores.shape[1] != 4:
        raise ValueError(
            f"{scores.shape} is incorrect shape for scores. "
            + "Expected (n_seqs, 4, seq_len)."
        )
    # save profile attributions
    np.savez_compressed(scores_path, scores * onehot_seqs)
    # save onehot seqs
    if onehot_seqs_path is not None:
        np.savez_compressed(onehot_seqs_path, onehot_seqs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "model_dir", type=str, help="Path to directory with model to score"
    )
    parser.add_argument(
        "fasta_path", type=str, help="Path to fasta file with sequences to score"
    )
    parser.add_argument("scores_path", type=str, help="where to save deepshap scores")
    parser.add_argument(
        "--ohe_seqs_path",
        type=str,
        default=None,
        help="where to save one-hot encoded sequences (optional, for TF-MoDISco)",
    )
    parser.add_argument("--n_shuffles", type=int, default=20)
    parser.add_argument("--rand", type=int, default=47)
    parser.add_argument(
        "--mode",
        type=str,
        help="Mode for attributions. Options are profile and counts",
        default="counts",
    )
    parser.add_argument("--gpu", type=int, default=None)
    parser.add_argument("--silence_tqdm", action="store_true")
    args = parser.parse_args()

    if args.gpu is not None:
        print(f"Using GPU {args.gpu}")
        device = f"cuda:{args.gpu}"
    else:
        device = "cpu"

    # Load model
    model_path = glob.glob(f"{args.model_dir}/*.torch")[0]
    loader = procapnet.ProCapNet()
    loader.load_state_dict(torch.load(model_path))
    if args.mode == "profile":
        # raise NotImplementedError("Profile mode not implemented")
        model = ProfileWrapper(loader)
    elif args.mode == "counts":
        model = CountWrapper(loader)
    else:
        raise ValueError(f"Mode {args.mode} not recognized")
    if args.gpu is not None:
        model = model.cuda()

    # Load data
    sequences = pyfastx.Fasta(args.fasta_path)
    ohe = torch.stack(
        [
            one_hot_encode(rec.seq.upper())
            for rec in tqdm.tqdm(
                sequences,
                desc="One-hot encoding",
                total=len(sequences),
                disable=args.silence_tqdm,
            )
        ],
    ).to(torch.float32)
    if args.gpu is not None:
        ohe = ohe.cuda()

    # Calculate attributions
    attributions = deep_lift_shap(
        model,
        ohe,
        n_shuffles=args.n_shuffles,
        device=device,
        batch_size=8,
        random_state=args.rand,
        verbose=not args.silence_tqdm,
    )

    # Save outputs
    save_deepshap_results(
        ohe.cpu().numpy(),
        attributions,
        args.scores_path,
        onehot_seqs_path=args.ohe_seqs_path,
    )


if __name__ == "__main__":
    main()
