## FTP 1000 genomes mitochondrial sequences
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz -P data/ThousandGenomes/
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz.tbi -P data/ThousandGenomes/

## Normalize indels, split multialleics, fill info
bcftools norm -f data/rCRS.fasta -m - \
data/ThousandGenomes/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz | \
bcftools view -V indels,mnps | \
bcftools norm -m + | \
bcftools +fill-tags -Oz -o data/ThousandGenomes/ALL.chrMT.vcf.gz

## Decompose
vt decompose data/ThousandGenomes/ALL.chrMT.vcf.gz | \
bcftools +fill-tags -Oz -o data/ThousandGenomes/ALL.chrMT_decomposed.vcf.gz

## Keep Only first alternate
./scripts/PYTHON/pickFirstAlt data/ThousandGenomes/ALL.chrMT_decomposed.vcf.gz | bgzip > data/ThousandGenomes/ALL.chrMT_decomposed_noMulti.vcf.gz
