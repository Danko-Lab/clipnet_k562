# PINTS: Peak Identifier for Nascent Transcripts Sequencing
![](https://img.shields.io/badge/platform-linux%20%7C%20osx-lightgrey.svg)
![](https://img.shields.io/badge/python-3.5%20%7C%203.6%20%7C%203.7-blue.svg)

## Prerequisite
Python packages
* matplotlib
* numpy
* pandas
* pysam
* pybedtools
* pyBigWig
* scipy
* statsmodels

Binary packages
* bedtools
* bgzip
* tabix


## Get started
PINTS can call peaks directly from BAM files. To call peaks from BAM files, 
you need to provide the tool a path to the bam file and what kind of experiment it was from.
If it's from a standard protocol, like PROcap, then you can set `--exp-type PROcap`.  
Other supported experiments including PROseq/PROcap/GROseq/GROcap/CoPRO. If the data was generated
by other methods, then you need to tell the tool where it can find ends of RNAs which you are interested in.
For example, `--exp-type R_5` tells the tool that:
   1. this is a single-end library; 
   2. the tool should look at 5' of reads. Other supported values are `R_3`, `R1_5`, `R1_3`, `R2_5`, `R2_3`.

If reads represent the reverse complement of original RNAs, like PROseq, then you need to use `--reverse-complement` 
(not necessary for standard protocols).

One example for calling peaks from BAM file:
```shell
./caller.py input.bam output_dir output_prefix --thread 16 --exp-type PROcap
```
Or you can call peaks from BigWig files:
```shell
./caller.py bigwig output_dir output_prefix --bw-pl path_to_pl.bw --bw-mn path_to_mn.bw --thread 16
```

## Outputs
* prefix+`_divergent_peaks.bed`: Divergent TREs;
* prefix+`_bidirectional_peaks.bed`: Bidirectional TREs (divergent + convergent);
* prefix+`_single_peaks.bed`: Unidirectional TREs, maybe lncRNAs transcribed from enhancers (e-lncRNAs) as suggested [here](http://www.nature.com/articles/s41576-019-0184-5).
* prefix+`_narrowPeak.bed`: All significant peaks in [narrowPeak](http://genome.ucsc.edu/FAQ/FAQformat.html#format12) format.

For divergent or bidirectional TREs, there will be 5 columns in the outputs:
1. Chromosome
2. Start site: 0-based
3. End site: 0-based
4. Confidence about the peak pair. Can be: 
    * `Stringent(qval)`, which means the two peaks on both forward and reverse strands are significant based-on their q-values; 
    * `Stringent(pval)`, which means one peak is significant according to q-value while the other one is significant according to p-value; 
    * `Relaxed`, which means only one peak is significant in the pair.
    * A combination of the three types above, because of overlap for nearby elements.
5. Peak ID

For single TREs, there will be  columns in the output:
1. Chromosome
2. Start
3. End
4. Peak ID
5. Q-value
6. Strand

## Parameters
### Required parameters
* `bam_file`: input bam file, if you want to use bigwig files, please leave this as bigwig
* `save_to`: save peaks to
* `file_prefix`: prefix to all intermediate files

### Conditionally required parameters
* `--reverse-complement`: Set this switch if reads in this library represent the reverse complement of nascent RNAs, like PROseq
* `--exp-type <experiment type>`: Type of experiment, acceptable values are: `CoPRO`/`GROcap`/`GROseq`/`PROcap`/`PROseq`, or if you know the position of RNA ends which you're interested on the reads, you can specify `R_5`, `R_3`, `R1_5`, `R1_3`, `R2_5` or `R2_3`;
* `--bw-pl <path to bigwig (pl)>`: Bigwig for plus strand. If you want to use bigwig instead of BAM, please set bam_file to bigwig;
* `--bw-mn <path to bigwig (mn)>`: Bigwig for minus strand. If you want to use bigwig instead of BAM, please set bam_file to bigwig;

### Optional parameters
* `--mapq-threshold <min mapq>`: Minimum mapping quality, by default: 30 or `None`;
* `--close-threshold <close distance>`: Distance threshold for two peaks (on opposite strands) to be merged, by default: 300;
* `--window-size <window size>`: Size for sliding windows, by default: 100;
* `--fdr-target <fdr>`: FDR target for multiple testing, by default: 0.1;
* `--adjust_method <fdr method>`: Method used for testing and adjustment of p-values. By default, Benjamini/Hochberg (fdr_bh), ohter options are available at [here](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html);
* `--chromosome-start-with <chromosome prefix>`: Only keep reads mapped to chromosomes with this prefix, if it's set to `None`, then all reads will be analyzed;
* `--thread <n thread>`: Max number of threads the tool can create;
* `--output-diagnostic-plot`: Save diagnostic plots (independent filtering and pval dist) to local folder

More parameters can be seen by running `./caller.py -h`.

## Other tools
* `element_boundary.py`: Extend peaks from summits.
* `visualizer.py`: Generate bigwig files for the inputs.
* `normalizer.py`: Normalize inputs.

## Tips
1. Be cautious to reads mapped to scaffolds instead of main chromosome (for example the notorious `chrUn_gl000220` in `hg19`, they maybe rRNA contamination)!

## Contact
Please submit an issue with any questions or if you experience any issues/bugs.
