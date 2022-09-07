#### Step1: Randomly select 50K variants from the whole genome as TF-occupied variants; Take 100th repeat as an example:
python Random_TF.py 100

#### Step2: Use the GWAS summary statistics to prioritize TF-occupied regulatory variants. We run the linear mixed model to test the significant of the TF:
Rscript TFselect.R 100

#### Step3: If the p-value of the TF<0.05, we include their corresponding 50K variants as the prioritized TF-occupied regulatory variants to train the elastic net:
For each chromsome from 1 to 22, we run below command:
Rscript sTF_TWAS.R "$chrom" 100

#### Step4: Combined the weights from each chromosome, and performed TWAS analysis using the same GWAS summary statistics:
Associations.cmd
