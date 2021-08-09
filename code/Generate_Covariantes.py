import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

# Load the dataset
filename='Breast_Mammary_Tissue.v8.normalized_expression.no_sex.female.bed'
df = pd.read_csv(filename,sep="\t")
df_gene=df.iloc[:,4:]
df_gene_T=df_gene.T
file_name_gene='GeneExpression_Breast_Female_ID.csv'
df_gene_T.to_csv(file_name_gene,index=False,header=False)

# Generate first five principal components for each sample
pca = PCA(n_components=5)
X_pca=pca.fit_transform(df_gene_T)
X_pca_df = pd.DataFrame(X_pca).add_prefix('PCA-')
X_pca_df.to_csv("Breast_Mammary_Tissue_nosex_female.Principal_Components.csv",index=False)

# Load age information for each sample and combine the dataa with principal components
file_age='pheno_id_age.txt'
df_age = pd.read_csv(file_age,sep="\t")
combine=[df_age,X_pca_df]
df_combine = pd.concat(combine,axis=1)
df_result=df_combine.iloc[:,1:]
file_covariates='Covariates_Age_PC.csv'
df_result.to_csv(file_covariates,index=False,header=False)
