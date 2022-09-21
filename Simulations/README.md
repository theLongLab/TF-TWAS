### Under the GWAS/ 
#### 1) Randomly generated the phenotype values (0 or 1) independent of the genotype:
Random_Pheno.tfam contains randomly generated phenotype values

#### 2) Conducted logistic regression analysis to generate the GWAS summary statistics using the phenotype values and the genotype data from the 1000 Genomes Project:
gwas_association.cmd

### Under the Sim_TF/ 
#### 1) Randomly assigned 50K TF-occupied variants to a value “1” and the remaining variants to a value “0” ; Take 100th repeat as an example:
python Random_TF.py 100

#### 2) We then used generalized mixed models to estimate an association between the Chi-squared values (Y) and TF binding status of genetic variants (see Equation 1 in the manuscript):
Rscript TFselect.R 100

#### 3) We prioritized a set of variants based on the association with cancer risk at P < 0.05. We included their corresponding 50K variants as the prioritized TF-occupied regulatory variants to train the elastic net:
For each chromsome from 1 to 22, we run below command:
Rscript sTF_TWAS.R "$chrom" 100

#### 4) Combined the weights from each chromosome, and performed TWAS analysis using the same GWAS summary statistics:
Associations.cmd
