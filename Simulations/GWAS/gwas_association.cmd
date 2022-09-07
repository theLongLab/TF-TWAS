plink="/path/to/plink"

$plink --noweb --tped ALL.GRCh38.overlap_GTEx.clean.tped --tfam Random_Pheno.tfam --make-bed --out ALL.GRCh38.overlap_GTEx.clean

$plink --bfile ALL.GRCh38.overlap_GTEx.clean --pca --out ALL.GRCh38.overlap_GTEx.clean

$plink -bfile ALL.GRCh38.overlap_GTEx.clean --allow-no-sex --covar ALL.GRCh38.overlap_GTEx.clean.eigenvec --covar-number 2-12 --out ALL.GRCh38.overlap_GTEx.clean_beta --logistic beta --ci 0.95 --adjust

sed -e 's/\s\+/,/g' ALL.GRCh38.overlap_GTEx.clean_beta.assoc.logistic > ALL.GRCh38.overlap_GTEx.clean_beta.assoc.logistic_cp
mv ALL.GRCh38.overlap_GTEx.clean_beta_cp ALL.GRCh38.overlap_GTEx.clean_beta.assoc.logistic
