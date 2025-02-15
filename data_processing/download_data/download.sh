
#!/bin/bash

outdir=../../data/raw
mkdir -p $outdir

# Download reference genome

echo "Downloading reference genome (hg38) ..."
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz \
    -O ${outdir}/hg38.fa.gz
gunzip ${outdir}/hg38.fa.gz

# Download rDNA reference sequence

rDNA_fa=${outdir}/rDNA_human_U13369.1.fa
wget https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi\?tool\=portal\&save\=file\&log$\=seqview\&db\=nuccore\&report\=fasta\&id\=555853\&conwithfeat\=on\&withparts\=on\&hide-cdd\=on -O - | \
    sed '$ d' > $rDNA_fa

# Then, combine the two fasta files

combo_fa="$outdir/hg38.withrDNA.fa"
cat ${outdir}/hg38.fa "$rDNA_fa" > "$combo_fa"

# Download chromosome sizes

wget https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/@@download/GRCh38_EBV.chrom.sizes.tsv -O - | \
    grep -v "_alt" > ${outdir}/hg38.chrom.sizes

# Create a chromosome sizes file for hg38 + rDNA

cp "$outdir/hg38.chrom.sizes" "$outdir/hg38.withrDNA.chrom.sizes"
rDNA_length=`grep -v ">" "$rDNA_fa" | tr -d '\n' | wc -c`
echo -e "U13369.1\t$rDNA_length" >> "$outdir/hg38.withrDNA.chrom.sizes"

# Download PRO-cap data

echo "Downloading PRO-cap data ..."
# pl strand
wget -O $outdir/ENCFF994CSC.bigWig https://www.encodeproject.org/files/ENCFF994CSC/@@download/ENCFF994CSC.bigWig
wget -O $outdir/ENCFF580QEK.bigWig https://www.encodeproject.org/files/ENCFF580QEK/@@download/ENCFF580QEK.bigWig
# mn strand
wget -O $outdir/ENCFF253HUA.bigWig https://www.encodeproject.org/files/ENCFF253HUA/@@download/ENCFF253HUA.bigWig
wget -O $outdir/ENCFF328OOU.bigWig https://www.encodeproject.org/files/ENCFF328OOU/@@download/ENCFF328OOU.bigWig

# Done

echo "Done downloading raw data to ${outdir}."

