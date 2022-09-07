"$chrom" #from 1 to 22
Rscript /work/long_lab/jingni/project/TF_TWAS/Simulation/Simulation_Revision/Sim_GWAS/Sim_TF/sTF_TWAS.R "$chrom" 100

python Combine_Sum.py 100
Rscript WriteSQL.R 100
gzip ./summary/GTEx_100_covariances.txt
/path/to/SPrediXcan.py \
--model_db_path ./summary/GTEx_Pheno0_100.db \
--covariance ./summary/GTEx_100_covariances.txt.gz \
--gwas_folder ./GWAS/ \
--gwas_file_pattern ".*gz" \
--snp_column SNP_ID \
--effect_allele_column A1 \
--non_effect_allele_column A2 \
--beta_column Beta \
--pvalue_column P \
--output_file ./output/twas_100.csv
