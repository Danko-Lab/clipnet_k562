"""
Calculate contribution scores using shap.DeepExplainer.
Also contains some wrapper and utility functions used in scripts that make use of shap.
"""

import argparse
import gc
import glob
import logging
import os

import numpy as np
import pyfastx
import shap
import tqdm
import utils

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "4"
logging.getLogger("tensorflow").setLevel(logging.FATAL)
import tensorflow as tf

# This will fix an error message for running tf.__version__==2.5
shap.explainers._deep.deep_tf.op_handlers["AddV2"] = (
    shap.explainers._deep.deep_tf.passthrough
)
tf.compat.v1.disable_v2_behavior()


def load_seqs(
    fasta_fp, return_twohot_explains=True, background_fp=None, n_subset=100, seed=None
):
    np.random.seed(seed)
    seqs_to_explain = pyfastx.Fasta(fasta_fp)
    background_seqs = (
        seqs_to_explain if background_fp is None else pyfastx.Fasta(background_fp)
    )
    reference = [
        background_seqs[i]
        for i in np.random.choice(
            np.array(range(len(background_seqs))),
            size=min(n_subset, len(background_seqs)),
            replace=False,
        )
    ]
    shuffled_reference = [
        utils.kshuffle(rec.seq, random_seed=seed)[0] for rec in reference
    ]
    twohot_background = np.array(
        [utils.TwoHotDNA(seq).twohot for seq in shuffled_reference]
    )
    if return_twohot_explains:
        seqs_to_explain = np.array(
            [utils.TwoHotDNA(seq).twohot for seq in seqs_to_explain]
        )
    return seqs_to_explain, twohot_background


def create_explainers(model_fps, twohot_background, silence=False):
    models = [
        tf.keras.models.load_model(fp, compile=False)
        for fp in tqdm.tqdm(model_fps, desc="Loading models")
    ]
    explainers = []
    for model in tqdm.tqdm(models, desc="Creating explainers", disable=silence):
        explainers.append(
            shap.DeepExplainer((model.input, model.output), twohot_background)
        )
    return explainers


def calculate_scores(
    explainers, seqs_to_explain, batch_size=256, check_additivity=True, silence=False
):
    hyp_explanations = {i: [] for i in range(len(explainers))}
    for i, explainer in enumerate(explainers):
        desc = "Calculating explanations"
        if len(explainers) > 1:
            desc += f" for model fold {i + 1}"
        for j in tqdm.tqdm(
            range(0, len(seqs_to_explain), batch_size), desc=desc, disable=silence
        ):
            shap_values = explainer.shap_values(
                seqs_to_explain[j : j + batch_size], check_additivity=check_additivity
            )
            hyp_explanations[i].append(shap_values)
            gc.collect()

    concat_explanations = [
        np.concatenate([exp[0] for exp in hyp_explanations[k]], axis=0)
        for k in hyp_explanations.keys()
    ]

    if len(explainers) > 1:
        mean_explanations = np.array(concat_explanations).mean(axis=0)
    else:
        mean_explanations = concat_explanations[0]
    explanations = mean_explanations * seqs_to_explain
    return explanations, mean_explanations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta_fp", type=str, help="Fasta file path.")
    parser.add_argument("score_fp", type=str, help="Where to write DeepSHAP scores.")
    parser.add_argument(
        "--ohe_seq_fp", type=str, default=None, help="Where to write onehot sequences."
    )
    parser.add_argument(
        "--model_dir",
        type=str,
        default="ensemble_models/",
        help="Directory to load models from",
    )
    parser.add_argument(
        "--model_fp",
        type=str,
        default=None,
        help="Model file path. Use to calculate for a specific model fold. \
            Selecting this option will override --model_dir.",
    )
    parser.add_argument(
        "--hyp_attr_fp",
        type=str,
        default=None,
        help="Where to write hypothetical attributions.",
    )
    parser.add_argument(
        "--gpu",
        type=int,
        default=0,
        help="Index of GPU to use (starting from 0). Use --gpu -1 to use CPU.",
    )
    parser.add_argument(
        "--n_subset",
        type=int,
        default=100,
        help="Maximum number of sequences to use as background. \
            Default is 100 to ensure reasonably fast compute on large datasets.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for selecting background sequences.",
    )
    parser.add_argument(
        "--silence",
        action="store_true",
        help="Disables progress bars and other non-essential print statements.",
    )
    parser.add_argument(
        "--skip_check_additivity",
        action="store_true",
        help="Disables check for additivity of shap results.",
    )
    args = parser.parse_args()

    # Check arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if args.model_fp is None and args.model_dir is None:
        raise ValueError("Must specify either --model_fp or --model_dir.")

    # Load sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    seqs_to_explain, twohot_background = load_seqs(
        args.fasta_fp, n_subset=args.n_subset, seed=args.seed
    )

    # Create explainers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Set VRAM usage to growth and enable GPU usage
    if args.gpu is not None and args.gpu != -1:
        gpus = tf.config.list_physical_devices("GPU")
        gpu = gpus[args.gpu]
        tf.config.experimental.set_memory_growth(gpu, True)
        tf.config.set_visible_devices(gpu, "GPU")
    else:
        tf.config.set_visible_devices([], "GPU")

    if args.model_fp is not None:
        model_fps = [args.model_fp]
    else:
        model_fps = list(glob.glob(os.path.join(args.model_dir, "*.h5")))
    explainers = create_explainers(
        model_fps, twohot_background, args.silence or len(model_fps) == 1
    )

    # Calculate scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    explanations, mean_explanations = calculate_scores(
        explainers,
        seqs_to_explain,
        batch_size=256,
        silence=args.silence,
        check_additivity=not args.skip_check_additivity,
    )

    # Save scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Save DeepSHAP scores
    np.savez_compressed(args.score_fp, explanations.swapaxes(1, 2))
    # Convert twohot to onehot and save
    if args.ohe_seq_fp is not None:
        np.savez_compressed(
            args.ohe_seq_fp, (seqs_to_explain / 2).astype(int).swapaxes(1, 2)
        )
    # Save hypothetical attributions
    if args.hyp_attr_fp is not None:
        np.savez_compressed(args.hyp_attr_fp, mean_explanations.swapaxes(1, 2))


if __name__ == "__main__":
    main()
