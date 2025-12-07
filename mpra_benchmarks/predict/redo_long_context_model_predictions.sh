bedtools slop -i media-3_oligos_snps_cleaned_holdouts.bed.gz -b 262144 -g /home2/ayh8/data/hg19.chrom.sizes | \
    awk '$3 - $2 == 524289 {print}' | \
    sort-bed - | \
    bgzip > media-3_oligos_snps_cleaned_holdouts_524288bp.bed.gz

# Get SNPs that are active in K562
zcat media-4-K562_allelic_mpra.tsv.gz | awk -F"\t" '$9 == 1{print $1}' | sort > media-4-K562_allelic_mpra_activeK562.txt

join -1 1 -2 4 \
    media-4-K562_allelic_mpra_activeK562.txt \
    <(gunzip -c media-3_oligos_snps_cleaned_holdouts_524288bp.bed.gz | sort -k 4) | \
    awk 'OFS="\t" {print $2, $3, $4, $1}' | \
    sort-bed - | \
    bgzip > media-3_oligos_snps_cleaned_holdouts_524288bp_activeK562.bed.gz