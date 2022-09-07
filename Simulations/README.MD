### Under the GWAS/ 
#### 1) Randomly generated the phenotype values (0 or 1) independent of the genotype:
Random_Pheno.tfam contains randomly generated phenotype values

#### 2) Computed the GWAS summary statistics using the genotype from the 1000 Human Genomes Project and the random traits (0 or 1):
gwas_association.cmd

### Under the Sim_TF/ 
#### 1) Randomly select 50K variants from the whole genome as TF-occupied variants; Take 100th repeat as an example:
python Random_TF.py 100

#### 2) Use the GWAS summary statistics to prioritize TF-occupied regulatory variants. We run the linear mixed model to test the significant of the TF:
Rscript TFselect.R 100

#### 3) If the p-value of the TF<0.05, we include their corresponding 50K variants as the prioritized TF-occupied regulatory variants to train the elastic net:
For each chromsome from 1 to 22, we run below command:
Rscript sTF_TWAS.R "$chrom" 100

#### 4) Combined the weights from each chromosome, and performed TWAS analysis using the same GWAS summary statistics:
Associations.cmd
