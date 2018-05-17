bcftools norm -f ~/Dropbox/src/MitoImpute/data/rCRS.fasta \
-m -any -Oz -o ~/Dropbox/src/MitoImpute/data/ThousandGenomes/ALL.chrMT.vcf.gz \
~/Dropbox/src/MitoImpute/data/ThousandGenomes/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz

bcftools +fill-tags -Oz -o ~/Dropbox/src/MitoImpute/data/ThousandGenomes/ALL.chrMT2.vcf.gz ~/Dropbox/src/MitoImpute/data/ThousandGenomes/ALL.chrMT.vcf.gz
