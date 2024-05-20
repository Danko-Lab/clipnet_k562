### Code here is modified from Jacob's Schreiber's
### implementation of BPNet, called BPNet-lite:
### https://github.com/jmschrei/bpnet-lite/


import numpy as np
import torch

torch.backends.cudnn.benchmark = True


class Model(torch.nn.Module):
    """A basic BPNet model with stranded profile and total count prediction.
    This is a reference implementation for BPNet. The model takes in
    one-hot encoded sequence, runs it through:
    (1) a single wide convolution operation
    THEN
    (2) a user-defined number of dilated residual convolutions
    THEN
    (3a) profile predictions done using a very wide convolution layer
    AND
    (3b) total count prediction done using an average pooling on the output
    from 2 followed by a dense layer.

    This implementation differs from the original BPNet implementation in
    two ways:
    (1) A single log softmax is applied across both strands such that
    the logsumexp of both strands together is 0. Put another way, the
    two strands are concatenated together, a log softmax is applied,
    and the MNLL loss is calculated on the concatenation.
    (2) The count prediction task is predicting the total counts across
    both strands. The counts are then distributed across strands according
    to the single log softmax from 1.

    Parameters
    ----------
    model_save_prefix: str
        filepath to save model and performance metrics to.
    n_filters: int, optional
        The number of filters to use per convolution. Default is 64.
    n_layers: int, optional
        The number of dilated residual layers to include in the model.
        Default is 8.
    n_outputs: int, optional
        The number of outputs from the model. Generally either 1 or 2
        depending on if the data is unstranded or stranded. Default is 2.
    alpha: float, optional
        The weight to put on the count loss.
    trimming: int or None, optional
        The amount to trim from both sides of the input window to get the
        output window. This value is removed from both sides, so the total
        number of positions removed is 2*trimming.
    """

    def __init__(
        self,
        model_save_path,
        n_filters=64,
        n_layers=8,
        n_outputs=2,
        alpha=1,
        trimming=None,
    ):
        super(Model, self).__init__()
        self.n_filters = n_filters
        self.n_layers = n_layers
        self.n_outputs = n_outputs
        self.alpha = alpha
        self.trimming = trimming or 2**n_layers
        self.model_save_path = model_save_path
        self.train_metrics = []

        self.iconv = torch.nn.Conv1d(4, n_filters, kernel_size=21, padding=10)

        self.rconvs = torch.nn.ModuleList(
            [
                torch.nn.Conv1d(
                    n_filters, n_filters, kernel_size=3, padding=2**i, dilation=2**i
                )
                for i in range(1, self.n_layers + 1)
            ]
        )

        self.deconv_kernel_size = 75  # will need in forward() to crop padding
        self.fconv = torch.nn.Conv1d(
            n_filters, n_outputs, kernel_size=self.deconv_kernel_size
        )

        self.relus = torch.nn.ModuleList(
            [torch.nn.ReLU() for _ in range(0, self.n_layers + 1)]
        )
        self.linear = torch.nn.Linear(n_filters, 1)

    def forward(self, X):
        start, end = self.trimming, X.shape[2] - self.trimming

        X = self.relus[0](self.iconv(X))
        for i in range(self.n_layers):
            X_conv = self.relus[i + 1](self.rconvs[i](X))
            X = torch.add(X, X_conv)

        X = X[
            :,
            :,
            start - self.deconv_kernel_size // 2 : end + self.deconv_kernel_size // 2,
        ]

        y_profile = self.fconv(X)

        X = torch.mean(X, axis=2)
        y_counts = self.linear(X).reshape(X.shape[0], 1)

        return y_profile, y_counts

    def predict(self, X, batch_size=64, logits=False):
        with torch.no_grad():
            starts = np.arange(0, X.shape[0], batch_size)
            ends = starts + batch_size

            y_profiles, y_counts = [], []
            for start, end in zip(starts, ends):
                X_batch = X[start:end]

                y_profiles_, y_counts_ = self(X_batch)
                if not logits:  # apply softmax
                    y_profiles_ = self.log_softmax(y_profiles_)
                y_profiles.append(y_profiles_.cpu().detach().numpy())
                y_counts.append(y_counts_.cpu().detach().numpy())

            y_profiles = np.concatenate(y_profiles)
            y_counts = np.concatenate(y_counts)
            return y_profiles, y_counts

    def log_softmax(self, y_profile):
        y_profile = y_profile.reshape(y_profile.shape[0], -1)
        y_profile = torch.nn.LogSoftmax(dim=-1)(y_profile)
        y_profile = y_profile.reshape(y_profile.shape[0], self.n_outputs, -1)
        return y_profile
