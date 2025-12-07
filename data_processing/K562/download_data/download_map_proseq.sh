# Use fasterq-dump to download the following SRA files:

# G1	K562	374,946,808             GSM1480327              SRR1554311,SRR1554312   (Core et al. 2014)
# G3	K562	57,520,888              GSM3452725	            SRR8137173              dREG
# G5	K562	71,452,942              GSM2361442+GSM2361443   SRR4454567,SRR4454568   (Vihervaara et al. 2017)
# G6	K562	26,972,822+27,046,373   GSM2545324+GSM2545325	SRR5364303,SRR5364304   (Dukler et al. 2017)

# Merge replicates and rename files

cat SRR1554311.fastq SRR1554312.fastq | pigz -p 8 > K562_proseq_G1.fastq.gz
cat SRR8137173.fastq | pigz -p 8 > K562_proseq_G3.fastq.gz
cat SRR4454567.fastq SRR4454568.fastq | pigz -p 8 > K562_proseq_G5.fastq.gz
cat SRR5364303.fastq SRR5364304.fastq | pigz -p 8 > K562_proseq_G6.fastq.gz

# Copy over reference files

gunzip -c /fs/cbsubscb17/storage/data/short_read_index/hg38/bwa-0.7.13-r1126/hg38.rRNA.fa.gz ../../
rsync --progress /fs/cbsubscb17/storage/data/short_read_index/hg38/hg38.rRNA.chrom.sizes  ../../
rsync --progress /fs/cbsubscb17/storage/data/short_read_index/hg38/bwa-0.7.13-r1126/hg38.rRNA.fa.*  ../../

# Run proseq2.0 to process PRO-cap data

proseq2.0 -i K562_proseq_G1.fastq.gz -o K562_proseq_G1

i=1 #3,5,6
mkdir -p K562_proseq_G${i}
bash ~/proseq2.0/proseq2.0.bsh \
    -SE \
    -i ../../hg38.rRNA.fa \
    -c ../../hg38.rRNA.chrom.sizes \
    -I K562_proseq_G${i} \
    -O K562_proseq_G${i} \
    -T K562_proseq_G${i}_tmp \
    -P \
    --thread=8

