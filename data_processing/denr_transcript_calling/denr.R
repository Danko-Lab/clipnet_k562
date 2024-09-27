# Adapted from the DENR (https://github.com/CshlSiepelLab/DENR) vignette (https://rpubs.com/gzbyzyx/DENR).


# Install and load the required packages
# BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(tidyr)
library(DENR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
hg19_txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
main_chromosomes <- paste0("chr", c(as.character(1:22), "X", "Y"))
txdb <- keepSeqlevels(hg19_txdb, main_chromosomes, pruning.mode = "coarse")

# Load the transcript data
gr_ds <- GenomicFeatures::transcripts(txdb, c("tx_name", "gene_id"))
gr_ds$gene_id <- as.character(gr_ds$gene_id)

# Create unique keys for each transcript
starts <- GenomicRanges::start(gr_ds)
ends <- GenomicRanges::end(gr_ds)
gr_ds$key <- paste0(
    GenomicRanges::seqnames(gr_ds), ":",
    starts, ":", 
    ends, ":",
    GenomicRanges::strand(gr_ds), ":",
    gr_ds$tx_name
)

# Set the bin size and bigwig file paths
bsize <- 250
bwp <- "/local/workdir/James/PauseEvolution/data/human_K562/K562_LC_1-2_QC_end_all-merge.plus.bw"
bwm <- "/local/workdir/James/PauseEvolution/data/human_K562/K562_LC_1-2_QC_end_all-merge.minus.bw"

# Call transcript abundance using uniform profile
tq <- transcript_quantifier(
    gr_ds,
    bin_size = bsize,
    transcript_name_column = "key",
    gene_name_column = "gene_id",
    mask_start_bins = c(1, 1),
    mask_end_bins = c(1, 1)
)
tq <- add_data(tq = tq, bigwig_plus = bwp, bigwig_minus = bwm)
tq <- fit(tq = tq)
head(transcript_abundance(tq))

# Call transcript abundance using shape profile
tsp <- transcript_shape_profile(
    transcripts = gr_ds,
    bigwig_plus = bwp,
    bigwig_minus = bwm,
    bin_size = bsize,
    linear_head_length = 500,
    linear_tail_length = 500,
    min_transcript_length = 2000
)
tq_shape <- apply_shape_profile(tq, tsp)
tq_shape <- fit(tq_shape)
ta_shape <- transcript_abundance(tq_shape)

bed <- separate_wider_delim(
    transcript_abundance(tq),
    #ta_shape,
    cols="transcript_name",
    delim=":",
    names=c("chrom", "start", "end", "strand", "tx_name")
)
write.table(
    bed,
    file="./denr_transcript_abundance.bed",
    quote=FALSE,
    row.names=FALSE,
    col.names=FALSE,
    sep="\t"
)