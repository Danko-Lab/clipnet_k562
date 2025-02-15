# PRO-seq and pausing index model notes

These describe unpublished work that is not well documented (sorry), but I'm keeping this here for convenience.

We fine tune the PRO-cap model to PRO-seq data:

```bash
python calculate_fold_params_proseq.py $DATADIR $OUTDIR
GPU=0
for i in {1..9}; do python ft_k562_proseq.py $i $GPU; done
```

And then the PRO-seq models to predict pausing index at promoters (also not documented here). First, calculate pausing indices and get promoter sequences using pipelines in `data_processing/pausing_index` and `data_processing/pausing_index_sequence`.

```bash
for i in {1..9}; do python ft_k562_pausing.py $i; done
```

For the gene body density models:

```bash
for i in {1..9}; do python ft_k562_initiation.py $i; done
```

```bash
clipnet predict
    -f ../../data/k562/initiation/k562_initiation_centered.fa.gz \
    -o ../../predictions/k562/k562_initiation_centered.npz \
    -m ../clipnet_k562/models/clipnet_k562_initiation/ \
    -n 1 -v
```

Benchmark initiation predictions:

```python
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr

expt = pd.read_csv("../../data/k562/initiation/k562_initiation_G156.bed", sep="\t", header=None)
pred = np.load("../../predictions/k562/k562_initiation_centered.npz")["arr_0"]

data = pd.DataFrame(
    {"chrom": expt.iloc[:, 0], "expt": expt.iloc[:, 4], "pred": pred.squeeze()}
)
test = data[data.chrom.isin([f"chr{c}" for c in [9, 13, 20, 21]])]
pearsonr(np.log1p(test.expt), np.log1p(test.pred))
# PearsonRResult(statistic=0.7002800253241194, pvalue=0.0)
spearmanr(test.expt, test.pred)
# SignificanceResult(statistic=0.6981647099246578, pvalue=0.0)
data.to_csv("../../predictions/k562/k562_initiation_G156_benchmark.csv.gz")
```

Calculate attributions for the pausing & initiation models:

```bash
time clipnet attribute \
    -f ../../data/k562/pausing_index/k562_pausing_index_centered.fa.gz \ 
    -o ../../attr/k562/k562_pausing_index_centered_deepshap.npz \ 
    -s ../../attr/k562/k562_pausing_index_centered_ohe.npz \
    -m ../models/clipnet_k562_pausing/ \
    -v -y

time clipnet attribute \
    -f ../../data/k562/initiation/k562_initiation_centered.fa.gz \
    -o ../../attr/k562/k562_initiation_centered_deepshap.npz \
    -s ../../attr/k562/k562_initiation_centered_ohe.npz \
    -m ../models/clipnet_k562_initiation/ \
    -v -y
```

Average DeepSHAP files, then calculate MODISCO:

```bash
cd ../../data/pausing_index/
time modisco motifs -s k562_pausing_index_centered_ohe.npz -a k562_pausing_index_centered_deepshap.npz -n 1000000 -l 50 -v -o k562_pausing_index_centered_modisco.h5
time modisco report -i k562_pausing_index_centered_modisco.h5 -o k562_pausing_index_centered_modisco/ -m /home2/ayh8/data/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
```
