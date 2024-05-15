outdir=../../data/raw/

mkdir -p $outdir

# Download reference genome

echo "Downloading reference genome (hg38) ..."
wget -O $outdir/hg38.fa.gz http://hg

# Download PRO-cap data

echo "Downloading PRO-cap data ..."
wget -O $outdir/PRO-cap_data.tar.gz http://hg

# Done

echo "Done downloading raw data to ${outdir}."