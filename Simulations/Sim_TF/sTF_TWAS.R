library(methods)
library(glmnet)
library(stringr)
library(data.table)
library(dplyr)

setwd("/path/to/work_dir/") 
#source("gtex_v7_nested_cv_elnet.R")
"%&%" <- function(a,b) paste(a,b, sep = "")

args = commandArgs(trailingOnly=TRUE)

chrom<-args[1]
index<-args[2]
set<-"50k"
prefix<-"GTEx_" %&% index
#### I. data file preparation #### 
snp_annot_file <- './SNP_Anno_"%&%set%&%"_" %&% index %&% '.txt'### need to loaded for top 50k/all_snp
gene_annot_file <- "/path/to/gencode.v26.GRCh38.genes.gtf"
genotype_file <- "/path/to/GTEx/Whole_Blood/GTEx_v8_866Indiv_wb.clean_chr"%&% chrom %&% ".num.csv" ## load csv file from GTEx, with the CHR,POS...
covariance_genotype_file <-"/path/to/GTEx/Whole_Blood/GTEx_v8_866Indiv_wb.clean_chr" %&% chrom %&% ".num.csv"
expression_file <- "/path/to/GTEx/Whole_Blood/Whole_Blood.v8.residuals.csv" ## gene epxression, have removed the PEER factors
covariates_file <- "GTEx_Cov_chr"%&% chrom %&% ".txt" ## To be notice: we used all SNP to calculate covariate
gwas_file <- "/TF_Binding_b38_"%&%index %&%".txt"

#### II. load data  ######
# snp_annot_file from GTEx-those variants that rank as top 50k
# Make different snp_anno_file for 50k,100k...all_snp, contain all the overlapping variants between GTEx and GWAS that column with
snp_annot <- read.table(snp_annot_file, header = T, stringsAsFactors = F)
snp_annot$chrpos<-snp_annot$SNP%&%"_"%&%snp_annot$ref%&%"_"%&%snp_annot$effect
gwassnp <- read.table(gwas_file,head=TRUE) ### need to visit for each chr by xingyi guo
gwassnp$chrpos<-gwassnp$SNP_ID%&%"_"%&%gwassnp$A1%&%"_"%&%gwassnp$A2
snp_annot<-snp_annot[snp_annot$chrpos %in% gwassnp$chrpos,]
print(dim(snp_annot))

get_gene_annotation <- function(gene_annot_file_name, chrom)
{
  gene_df <- read.table(gene_annot_file,header=F,stringsAsFactors =F,sep="\t",fill=T)
  gene_df1 <- filter(gene_df,V3 %in% "gene")
  geneid <- str_extract(gene_df1[,9], "ENSG\\d+.\\d+")
  genename <- gsub("gene_name (\\S+);","\\1",str_extract(gene_df1[,9], "gene_name (\\S+);"), perl=T)
  gene_used <- as.data.frame(cbind(geneid,genename,gene_df1[,c(1,4,5,3)]))
  colnames(gene_used) <- c("geneid","genename","chr","start","end","anno")
  gtf_used <- filter(gene_used,gene_used[,3] %in% ('chr' %&% chrom))
  gtf_used
}

get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  #expr_df <- as.data.frame(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = 1))
  expr_df <- as.data.frame(read.table(file=gene_expression_file_name,header=FALSE,sep=",",check.names = FALSE))
  col_name=expr_df$V4
  expr_df<- expr_df[,c(-1,-2,-3,-4)]
  expr_df_T = as.data.frame(t(expr_df))
  colnames(expr_df_T)=col_name
  rownames(expr_df_T)=colnames(expr_df)
  expr_df_T <- expr_df_T %>% select(one_of(intersect(gene_annot$geneid, colnames(expr_df_T))))
  expr_df_T
}

get_filtered_genotype <- function(genotype_file_name) {
  gt_df <- fread(file=genotype_file_name,header=T,sep=",")
  gt_df <- as.data.frame(gt_df)
  snp_name<-gt_df$CHR%&%"_"%&%gt_df$LOC
  gt_df_T<- as.data.frame(t(as.matrix(gt_df[, c(-1,-2)])))
  colnames(gt_df_T) <- snp_name 
  gt_df_T
}

get_covariance_genotype<-function(covariance_genotype_file) {
  gt_df <- fread(file=covariance_genotype_file,header=T,sep=",")
  gt_df <- as.data.frame(gt_df)
  snp_name<-gt_df$CHR%&%"_"%&%gt_df$LOC
  gt_df_T<- as.data.frame(t(as.matrix(gt_df[, c(-1,-2)])))
  colnames(gt_df_T) <- snp_name
  gt_df_T
}

get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$geneid == gene),]
  c(row$start, row$end)
}

get_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {
  snp_info <- snp_annot %>% filter((pos >= (coords[1] - cis_window)  & (pos <= (coords[2] + cis_window))))
  cis_gt <- gt_df[, colnames(gt_df) %in% snp_info$SNP]
  cis_gt
}

do_elastic_net <- function(cis_gt, expr_adj, n_folds, cv_fold_ids, n_times, alpha) {
  #tryCatch({
  cis_gt <- as.matrix(cis_gt)
  fit <- cv.glmnet(cis_gt, expr_adj, nfolds = n_folds, alpha = alpha, keep = TRUE, type.measure='mse', foldid = cv_fold_ids[,1], parallel = FALSE)
  lambda_seq <- fit$lambda
  cvms <- matrix(nrow=length(lambda_seq), ncol=n_times)
  fits <- list()
  fits[[1]] <- fit
  cvms <- matrix(nrow = 100, ncol = n_times)
  cvms[1:length(fit$cvm),1] <- fit$cvm
  for (i in 2:(n_times)) {
    fit <- cv.glmnet(cis_gt, expr_adj, lambda = lambda_seq, nfolds = n_folds, alpha = alpha, keep = FALSE, foldid = cv_fold_ids[,i], parallel = FALSE)
    fits[[i]] <- fit
    cvms[1:length(fit$cvm),i] <- fit$cvm
  }
  avg_cvm <- rowMeans(cvms)
  best_lam_ind <- which.min(avg_cvm)
  best_lambda <- lambda_seq[best_lam_ind]
  out <- list(cv_fit = fits[[1]], min_avg_cvm = min(avg_cvm, na.rm = T), best_lam_ind = best_lam_ind, best_lambda = best_lambda)
  #  out
  #},error = function(cond) {
  #  message('Error')
  #  message(geterrmessage())
  #  out <- list()
  #  out
  #}
  #)
  out
}

evaluate_performance <- function(cis_gt, expr_adj, fit, best_lam_ind, best_lambda, cv_fold_ids, n_folds) {
  n_nonzero <- fit$nzero[best_lam_ind]
  if (n_nonzero > 0) {
    R2 <- rep(0, n_folds)
    for (j in (1:n_folds)) {
      fold_idxs <- which(cv_fold_ids[,1] == j)
      tss <- sum(expr_adj[fold_idxs]**2)
      rss <- sum((expr_adj[fold_idxs] - fit$fit.preval[fold_idxs, best_lam_ind])**2)
      R2[j] <- 1 - (rss/tss)
    }
    best_fit <- fit$glmnet.fit
    expr_adj_pred <- predict(best_fit, as.matrix(cis_gt), s = best_lambda)
    tss <- sum(expr_adj**2)
    rss <- sum((expr_adj - expr_adj_pred)**2)
    
    n_samp <- length(expr_adj)
    #f_stat <- ((tss - rss) / n_nonzero) / (rss / (n_samp - ncol(cis_gt) - 1))
    #p_val <- pf(f_stat, n_samp - 1, n_samp - ncol(cis_gt) - 1, lower.tail = FALSE)
    weights <- best_fit$beta[which(best_fit$beta[,best_lam_ind] != 0), best_lam_ind]
    weighted_snps <- names(best_fit$beta[,best_lam_ind])[which(best_fit$beta[,best_lam_ind] != 0)]
    R2_mean <- mean(R2)
    R2_sd <- sd(R2)
    inR2 <- 1 - (rss/tss)
    # Old way
    pred_perf <- summary(lm(expr_adj ~ fit$fit.preval[,best_lam_ind]))
    pred_perf_rsq <- pred_perf$r.squared
    
    pred_perf_pval <- pred_perf$coef[2,4]
    #one_sided_pval <- cor.test(expr_adj, fit$fit.preval[,best_lam_ind], alternative = 'greater')$p.value
    out <- list(weights = weights, n_weights = n_nonzero, weighted_snps = weighted_snps, R2_mean = R2_mean, R2_sd = R2_sd,
                inR2 = inR2, pred_perf_rsq = pred_perf_rsq, pred_perf_pval = pred_perf_pval)
  } else {
    out <- list(weights = NA, n_weights = n_nonzero, weighted_snps = NA, R2_mean = NA, R2_sd = NA,
                inR2 = NA, pred_perf_rsq = NA, pred_perf_pval = NA)
  }
  out
}

do_covariance <- function(gene_id, gt_cov_df, rsids, varIDs, out_file) {
  model_gt <- gt_cov_df[,varIDs, drop=FALSE]
  geno_cov <- cov(model_gt)
  #print(dim(geno_cov))
  cov_df <- data.frame(gene=character(),rsid1=character(),rsid2=character(), covariance=double())
  for (i in 1:length(rsids)) {
    for (j in i:length(rsids)) {
      #print(c(i, j))
      cov_df <- tryCatch(rbind(cov_df, data.frame(gene=gene_id,rsid1=rsids[i], rsid2=rsids[j], covariance=geno_cov[i,j])),
                         error = function(cond) browser())
    }
  }
  write.table(cov_df, file = out_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
}


### main function #####

n_times=3
n_k_folds=10
cis_window=1000000
alpha=0.5
# Read in data----
gene_annot <- get_gene_annotation(gene_annot_file, chrom)
print(dim(gene_annot))
expr_df <- get_gene_expression(expression_file, gene_annot)
print("expr_df")
print(dim(expr_df))
samples <- rownames(expr_df)
print("samples")
print(samples[1:10])
n_samples <- length(samples)
print(n_samples)
genes <- colnames(expr_df)
print(genes[1:10])
n_genes <- length(expr_df)
print(n_genes)
gt_df <- get_filtered_genotype(genotype_file)
print(dim(gt_df))
gt_cov_df<-get_covariance_genotype(covariance_genotype_file)
print("gt_cov_df")
print(dim(gt_cov_df))

seed <- NA
seed <- ifelse(is.na(seed), sample(1:1000000, 1), seed)
set.seed(seed)
cv_fold_ids <- matrix(nrow = n_samples, ncol = n_times)
for (j in 1:n_times)
  cv_fold_ids[,j] <- sample(1:n_k_folds, n_samples, replace = TRUE)
#  covariates_df <- get_covariates(covariates_file, samples)

#   Set seed and cross-validation fold ids----
#  set.seed(10001)

# Prepare output data----
model_summary_file <- './summary/' %&% prefix %&% '_chr' %&% chrom %&% '_model_summaries.txt'
model_summary_cols <- c('gene_id', 'gene_name', 'alpha', 'cv_mse', 'lambda_iteration', 'lambda_min', 'n_snps_in_model',
                        'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2', 'pred_perf_R2', 'pred_perf_pval')
write(model_summary_cols, file = model_summary_file, ncol = 12, sep = '\t')

weights_file <- './summary/' %&% prefix %&% '_chr' %&% chrom %&% '_weights.txt'
weights_col <- c('gene_id', 'rsid', 'varID', 'ref', 'alt', 'beta')
write(weights_col, file = weights_file, ncol = 6, sep = '\t')

# tiss_chr_summ_f <- '../summary/' %&% prefix %&% '_chr' %&% chrom %&% '_tiss_chr_summary.txt'
# tiss_chr_summ_col <- c('n_samples', 'chrom', 'cv_seed', 'n_genes')
# tiss_chr_summ <- data.frame(n_samples, chrom, seed, n_genes)
# colnames(tiss_chr_summ) <- tiss_chr_summ_col
# write.table(tiss_chr_summ, file = tiss_chr_summ_f, quote = FALSE, row.names = FALSE, sep = '\t')

covariance_file <- './summary/' %&% prefix %&% '_chr' %&% chrom %&% '_covariances.txt'
covariance_col <- c('GENE', 'RSID1', 'RSID2', 'VALUE')
write(covariance_col, file = covariance_file, ncol = 4, sep = ' ')

for (i in 1:n_genes) {
  #cat(i, "/", n_genes, "\n")
  gene <- genes[i]
  #print(gene)
  gene_name <- gene_annot$genename[gene_annot$geneid == gene]
  model_summary <- c(gene, gene_name, alpha, NA, NA, NA, 0, NA, NA, NA, NA, NA)
  coords <- get_gene_coords(gene_annot, gene)
  cis_gt <- get_cis_genotype(gt_df, snp_annot, coords, cis_window)
  print(dim(cis_gt))
  if(length(dim(cis_gt))==2){
        if (ncol(cis_gt) >= 2) {
                adj_expression <- expr_df[,i]
                #adj_expression <- adjust_for_covariates(expression_vec, covariates_df)
                elnet_out <- do_elastic_net(cis_gt, adj_expression, n_k_folds, cv_fold_ids, n_times, alpha)
                #print(elnet_out)
                if (length(elnet_out) > 0) {
                        eval <- evaluate_performance(cis_gt, adj_expression, elnet_out$cv_fit, elnet_out$best_lam_ind, elnet_out$best_lambda, cv_fold_ids, n_k_folds)
                        model_summary <- c(gene, as.character(gene_name), alpha, elnet_out$min_avg_cvm, elnet_out$best_lam_ind,
                        elnet_out$best_lambda, eval$n_weights, eval$R2_mean, eval$R2_sd, eval$inR2,
                        eval$pred_perf_rsq, eval$pred_perf_pval)
                        print("elnet_out")
                        #return(list(eval = eval, model_summary = model_summary, elnet_out = elnet_out, snp_annot = snp_annot))
                        if (eval$n_weights > 0) {
                                print(eval$weighted_snps)
                                print(eval$weights)
                                #print(colnames(snp_annot))
                                weighted_snps_info <- snp_annot %>% filter(SNP %in% eval$weighted_snps) %>% select(SNP,chrpos,ref,effect)
                                print(weighted_snps_info)
                                if (nrow(weighted_snps_info) != 0){
                                        weighted_snps_info$gene <- gene
                                        weighted_snps_info <- weighted_snps_info %>% merge(data.frame(weights = eval$weights, SNP=eval$weighted_snps), by = 'SNP') %>%
                                        select(gene, SNP, chrpos, ref, effect, weights)
                                        write.table(weighted_snps_info, file = weights_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
                                        do_covariance(gene, gt_cov_df, weighted_snps_info$chrpos, weighted_snps_info$SNP, covariance_file)
                                }else if (nrow(weighted_snps_info) == 0){
                                        model_summary <- c(gene, gene_name, alpha, NA, NA, NA, 0, NA, NA, NA, NA, NA)
                                }
                        }
                }
        }
  }
  write(model_summary, file = model_summary_file, append = TRUE, ncol = 12, sep = '\t')
}

