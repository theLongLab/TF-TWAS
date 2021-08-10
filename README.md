## Users’ Manual of TF-TWAS

## Overview
This directory contains how we proccessed Elastic-net training and Generalized Berk-Jones (GBJ) cross-set association test that used in the paper titled "Integrating hierarchical transcription factors (TF)-occupied elements with transcriptome-wide association analysis identifies 121 putative susceptibility genes in breast cancer".

## Methods
### Elastic-Net training: 
#### Prepare input data: 
**1)	Gene expression file:** 

The original gene expression data “Breast_Mammary_Tissue.v8.normalized_expression.bed” was downloaded from GTEx portal(https://gtexportal.org/home/datasets). We only included 151 female samples in our analysis and generated a new file named “Breast_Mammary_Tissue.v8.normalized_expression.no_sex.female.bed”. This below command will generate a file named “Covariates_Age_PC.csv”, which contains covariates matrix of the top five PCs and age for each sample (N × C matrix, where N is the number of samples and C is the number of covariates).

`python ./code/Generate_Covariantes.py`

Then, we further performed a probabilistic estimation of expression residuals (PEER) analysis to adjust for top principal components (PCs), age, and other potential confounding factors (PEERs)[1] for downstream prediction model building. There is a description of how to download and use the PEER tool here: https://github.com/PMBio/peer/wiki/Tutorial. The command that we used is shown as below: 

`Rscript ./code/Peer_Script.R`

According to GTEx protocol, if the number of samples is between 150 and 250, 30 PEER factors should be used. For our study, the number of samples is 151, so we used 30 PEER factors. This command will generate a residual file named “GeneExpression_Breast_Female_AfterRM_Residuals.csv”, and from this residual file, we generated the final gene expression data file named “Breast_Mammary_Tissue.v8.normalized_expression.no_sex.female.rm_covariates.bed” as the input for our downstream predictive model. 

**2)	genotype file:**  
*.breast.female.no_sex.csv for each gene ID. This file only includes variants within +/- 1M of the target gene location. The format of the genotype file looks like below: 

SNP,CHR,LOC,GTEX-1117F,GTEX-1122O,GTEX-11EM3,GTEX-11EMC,GTEX-11GSP,GTEX-11I78, … …
chr1_933303_C_T_b38,1,933303,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, … …
chr1_933411_C_T_b38,1,933411,1,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, … …
chr1_933653_C_T_b38,1,933653,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0, … …
chr1_933741_T_TG_b38,1,933741,0,2,2,1,2,1,2,1,0,2,1,2,0,1,2,0,2,1,2,1,2,2,1,1,1,2,1,2,2,1,2,1,1 … …
… …


#### Training Elastic Net model:
Once every input file is ready, we processed one gene at a time by executing this code:	

`Rscript ./code/WeightEN_Stratify_GTexv7_New.R  target_gene_id`


### Generalized Berk-Jones (GBJ) cross-set association test:
In order to combine single-set test results, we adapted existing codes for integrating correlated and sparse signals in multiple tissues provided by UTMOST[2] (https://github.com/Joker-Jerome/UTMOST). First, a covariance matrix was built based on the single-set test results and LD structures are formed, then the GBJ test is carried out. The codes for both steps are shown as below:

**Step 1: Calculate the joint test covariance**

`python /path/to/joint_covariance.py --weight_db weight_db/ --input_folder input_folder/ --covariance_output covariance/`

**Command parameters:**


joint_covariance.py. It's downloaded from https://github.com/Joker-Jerome/UTMOST

weight_db. Name of weight db in data folder.

input_folder. Name of folder containing dosage data.

covariance_output. Path where covariance results will be saved to.


**Step 2: Combine gene-trait associations in variant set models by joint GBJ test**

`python /path/to/joint_GBJ_test.py --weight_db weight_db/ --output_dir result_gbj/ --cov_dir covariance/ --input_folder asso_test/ --gene_info intermediate/gene_info.txt --output_name Func_GBJ_Joint --start_gene_index 1 --end_gene_index 12700\`

**Command parameters:**

joint_GBJ_test.py. It's downloaded from https://github.com/Joker-Jerome/UTMOST

weight_db. Name of weight db in data folder (imputation models).

input_folder. Name of folder containing single-test association results.

cov_dir. Path where covariance results are (covariance matrix for gene-level test statistics across each set).

output_dir. Path where results will be saved to.

gene_info. File containing the all the genes tested.

start_gene_index. Index of the starting gene in intermediate/gene_info.txt.

end_gene_index. Index of the ending gene in intermediate/gene_info.txt. 

### References: 
1. Stegle, O., Parts, L., Piipari, M., Winn, J., and Durbin, R. (2012). Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses. Nat Protoc 7, 500-507.
2. Hu, Y., Li, M., Lu, Q., Weng, H., Wang, J., Zekavat, S.M., Yu, Z., Li, B., Gu, J., Muchnik, S., et al. (2019). A statistical framework for cross-tissue transcriptome-wide association analysis. Nat Genet 51, 568-576.

### Contacts
  Jingni He: jingni.he1@ucalgary.ca<br>
  Quan Long: quan.long@ucalgary.ca<br>
  Xingyi Guo: xingyi.guo@vumc.org<br>
