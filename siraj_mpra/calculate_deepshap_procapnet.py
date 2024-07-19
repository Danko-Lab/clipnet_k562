import argparse
import glob
import sys

import numpy as np
import torch
from tangermeme.deep_lift_shap import deep_lift_shap

sys.path.append("../")
import utils


class ProfileWrapper(torch.nn.Module):
    def __init__(self, model):
        super(ProfileWrapper, self).__init__()
        self.model = model

    def forward(self, X):
        return self.model(X)[0]


class CountsWrapper(torch.nn.Module):
    def __init__(self, model):
        super(CountsWrapper, self).__init__()
        self.model = model

    def forward(self, X):
        return self.model(X)[1]


def get_attributions(sequences, model, n_shuffles=5, device="cpu", rand=47):
    if len(sequences.shape) != 3 or sequences.shape[1] != 4:
        raise ValueError(
            f"{sequences.shape} is incorrect shape for one-hot encoded sequences. "
            + "Expected (n_seqs, 4, seq_len)."
        )

    attrs_fwd = deep_lift_shap(
        model,
        sequences,
        n_shuffles=n_shuffles,
        device=device,
        batch_size=8,
        random_state=rand,
    )
    # attrs_rev = deep_lift_shap(model, torch.flip(sequences, [1, 2]), n_shuffles=n_shuffles, device=device, batch_size=8, random_state=rand)
    # attrs_rev = torch.flip(attrs_rev, [1, 2])
    # Append mean of attributions for fwd and rev strands
    # attrs_batch = np.array([attrs_fwd.numpy(), attrs_rev.numpy()])
    # attrs.append(attrs_batch.mean(axis=0))
    return attrs_fwd


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
    parser.add_argument("--n_shuffles", type=int, default=5)
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
    import torch

    # Load model
    model_path = glob.glob(f"{args.model_dir}/*.torch")[0]
    if args.mode == "profile":
        raise NotImplementedError("Profile mode not implemented")
        loader = torch.load(model_path)
        loader.load_state_dict(torch.load(model_path))
        model = ProfileWrapper(loader)
    elif args.mode == "counts":
        loader = torch.load(model_path)
        loader.load_state_dict(torch.load(model_path))
        model = CountsWrapper(loader)
    else:
        raise ValueError(f"Mode {args.mode} not recognized")
    if args.gpu is not None:
        model = model.cuda()

    # Load data
    sequences = (
        torch.tensor(
            utils.get_twohot_fasta_sequences(args.fasta_path), dtype=torch.float32
        ).swapaxes(1, 2)
        / 2
    )

    # Calculate attributions
    attributions = get_attributions(
        sequences,
        model,
        n_shuffles=args.n_shuffles,
        device=device,
        rand=args.rand,
    )

    # Save outputs
    save_deepshap_results(
        sequences, attributions, args.scores_path, onehot_seqs_path=args.ohe_seqs_path
    )


if __name__ == "__main__":
    main()
