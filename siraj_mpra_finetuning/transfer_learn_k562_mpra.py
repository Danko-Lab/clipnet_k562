import logging
import os
import sys
from pathlib import Path

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "4"
logging.getLogger("tensorflow").setLevel(logging.FATAL)
import tensorflow as tf
from tensorflow.keras.callbacks import CSVLogger, LearningRateScheduler
from tqdm.keras import TqdmCallback

sys.path.append("../../clipnet/")
import clipnet
import mpra_gen
import rnn_v10

fold = int(sys.argv[1])
gpu = int(sys.argv[2])


# Specify GPU usage
nn = clipnet.CLIPNET(n_gpus=1, use_specific_gpu=gpu)


def corr(x, y, pseudocount=1e-6):
    """
    Computes Pearson's r between x and y. Pseudocount ensures non-zero denominator.
    """
    mx = tf.math.reduce_mean(x)
    my = tf.math.reduce_mean(y)
    xm, ym = x - mx, y - my
    num = tf.math.reduce_mean(tf.multiply(xm, ym))
    den = tf.math.reduce_std(xm) * tf.math.reduce_std(ym) + pseudocount
    r = tf.math.maximum(tf.math.minimum(num / den, 1), -1)
    return r


def warmup_lr(epoch, lr):
    """
    Learning rate warmup schedule.
    """
    print(f"LEARNING RATE = {lr}")
    if epoch < 1:
        return lr / 100
    elif epoch == 1:
        return lr * 10
    else:
        return lr


# Create data loaders for training and validation
data_fp = "../../data/k562_mpra/processed_k562_mpra_data_clipnet_ft.csv.gz"

train_folds = [i for i in range(10) if i not in [fold, fold % 9 + 1, 0]]
val_folds = [fold % 9 + 1]
print(f"Training on folds {train_folds} and validating on fold {val_folds}.")

train_args = [data_fp, train_folds, rnn_v10.batch_size]
val_args = [data_fp, val_folds, rnn_v10.batch_size]
train_gen = mpra_gen.MPRAGen(*train_args)
val_gen = mpra_gen.MPRAGen(*val_args)

# Load the reference and alternative models
outdir = Path(f"../models/clipnet_k562_mpra/f{fold}/")

ref_model = tf.keras.models.load_model(
    f"../models/clipnet_k562/fold_{fold}.h5", compile=False
)
for layer in ref_model.layers:
    layer._name = layer.name + str("_ref")

alt_model = tf.keras.models.load_model(
    f"../models/clipnet_k562/fold_{fold}.h5", compile=False
)
for layer in alt_model.layers:
    layer._name = layer.name + str("_alt")

# Create a new model that outputs the log2 fold change
output = ref_model.output[1] - alt_model.output[1]
mpra_net = tf.keras.Model(inputs=[ref_model.input, alt_model.input], outputs=output)
tf.keras.backend.clear_session()
mpra_net.compile(
    optimizer=rnn_v10.optimizer(**rnn_v10.opt_hyperparameters),
    loss="mse",
    metrics=["mae", corr],
)
for layer in mpra_net.layers:
    if isinstance(layer, tf.keras.layers.BatchNormalization):
        layer.trainable = False
    else:
        layer.trainable = True

# Compile the model
mpra_net_filepath = str(outdir.joinpath("mpra_net{epoch:02d}.tf"))
chkpt = tf.keras.callbacks.ModelCheckpoint(
    mpra_net_filepath, verbose=0, save_best_only=True
)
early_stopping = tf.keras.callbacks.EarlyStopping(verbose=1, patience=10)
tqdm_callback = TqdmCallback(verbose=1, bar_format="{l_bar}{bar:10}{r_bar}{bar:-10b}")
csv_logger = CSVLogger(
    filename=outdir.joinpath("mpra_net.log"), separator=",", append=True
)

# Fit the model
fit_model = mpra_net.fit(
    x=train_gen,
    validation_data=val_gen,
    epochs=rnn_v10.epochs,
    # batch_size=rnn_v10.batch_size,
    steps_per_epoch=train_gen.steps_per_epoch,
    validation_steps=val_gen.steps_per_epoch,
    verbose=0,
    callbacks=[
        chkpt,
        early_stopping,
        tqdm_callback,
        csv_logger,
        LearningRateScheduler(warmup_lr),
    ],
)
