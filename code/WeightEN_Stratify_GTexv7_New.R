library(glmnet)
"%&%" <- function(a,b) paste(a,b, sep = "")

args = commandArgs(trailingOnly=TRUE)

target_gene_id<-args[1]

working_dir=paste('/path/to/', target_gene_id,"/",sep="")
if(!dir.exists(working_dir)){
  dir.create(working_dir)
}
output_dir= paste(working_dir,"output/",sep="")
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

data_dir <-paste('/path/to/',target_gene_id ,'/',sep="")
inCsvfile <- paste(data_dir,target_gene_id,".breast.female.no_sex.csv",sep="")
df_geno <- read.table(file=inCsvfile,header=TRUE,sep=",",check.names = FALSE)
geno_indicator<- df_geno$SNP

df_geno_T<- as.data.frame(t(as.matrix(df_geno[, c(-1,-2,-3)])))

inbedfile<-'/path/to/Breast_Mammary_Tissue.v8.normalized_expression.no_sex.female.rm_covariates.bed'
df_bed <- read.table(file=inbedfile,header=FALSE,sep="\t",check.names = FALSE)
df_target_gene<- df_bed[df_bed$V4==target_gene_id,c(-1,-2,-3,-4)]
df_target_gene_T = as.data.frame(t(df_target_gene))
colnames(df_target_gene_T)="Gene_Exp"
rownames(df_target_gene_T)=as.vector(rownames(df_geno_T))
adj_expression<-df_target_gene_T$Gene_Exp

calc_R2 <- function(y, y_pred) {
  tss <- sum(y**2)
  rss <- sum((y - y_pred)**2)
  1 - rss/tss
}

generate_fold_ids <- function(n_samples, n_folds=10) {
  n <- ceiling(n_samples / n_folds)
  fold_ids <- rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}

nested_cv_elastic_net_perf <- function(x, y, n_samples, n_train_test_folds, n_k_folds, alpha) {
  R2_folds <- rep(0, n_train_test_folds)
  corr_folds <- rep(0, n_train_test_folds)
  zscore_folds <- rep(0, n_train_test_folds)
  pval_folds <- rep(0, n_train_test_folds)
  best_lam_list <- rep(0, n_train_test_folds)
  # Outer-loop split into training and test set.
  train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
  fits<-list()
  for (test_fold in 1:n_train_test_folds) {
    train_idxs <- which(train_test_fold_ids != test_fold)
    test_idxs <- which(train_test_fold_ids == test_fold)
    x_train <- x[train_idxs, ]
    y_train <- y[train_idxs]
    x_test <- x[test_idxs, ]
    y_test <- y[test_idxs]
    # Inner-loop - split up training set for cross-validation to choose lambda.
    cv_fold_ids <- generate_fold_ids(length(y_train), n_k_folds)
    y_pred <- tryCatch({
      # Fit model with training data.
      fit <- cv.glmnet(x_train, y_train, nfolds = n_k_folds, alpha = alpha,type.measure='mse', foldid = cv_fold_ids)
      fits[[test_fold]]<-fit
      best_lam_list[test_fold] <- which.min(fit$cvm)
      # Predict test data using model that had minimal mean-squared error in cross validation.
      predict(fit, x_test, s = 'lambda.min')},
      # if the elastic-net model did not converge, predict the mean of the y_train (same as all non-intercept coef=0)
      error = function(cond) rep(mean(y_train), length(y_test)))
    R2_folds[test_fold] <- calc_R2(y_test, y_pred)
    # Get p-value for correlation test between predicted y and actual y.
    # If there was no model, y_pred will have var=0, so cor.test will yield NA.
    # In that case, give a random number from uniform distribution, which is what would
    # usually happen under the null.
    corr_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor(y_pred, y_test), 0)
    zscore_folds[test_fold] <- atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3) # Fisher transformation
    pval_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
  }
  R2_avg <- mean(R2_folds)
  best_fit=fits[[which(R2_folds==max(R2_folds))]]
  best_lamba_select=best_lam_list[which(R2_folds==max(R2_folds))]
  R2_sd <- sd(R2_folds)
  rho_avg <- mean(corr_folds)
  rho_se <- sd(corr_folds)
  rho_avg_squared <- rho_avg**2
  
  # Stouffer's method for combining z scores.
  zscore_est <- sum(zscore_folds) / sqrt(n_train_test_folds)
  zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
  # Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
  list(fit=best_fit,best_lambda=best_lamba_select,R2_avg=R2_avg, R2_sd=R2_sd, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval)
}

### main function #####
n_folds=10
n_times=3
alpha=0.5
cis_gt=as.matrix(df_geno_T)

adj_expression=as.vector(adj_expression)
n_samples=length(adj_expression)
n_folds=10
n_train_test_folds=5

seed <- ifelse(is.na(seed), sample(1:1000000, 1), seed)
set.seed(seed)

# Prepare output data----
model_summary_file <- output_dir %&% target_gene_id %&% '_model_summary_GTExv7.txt'
model_summary_cols <- c('gene_id', 'alpha', 'n_snps_in_window', 'n_snps_in_model', 'n_snps_in_train_model', 'lambda_min_mse',
                        'test_R2_avg', 'test_R2_sd', 'cv_R2_avg', 'cv_R2_sd', 'training_R2', 'all_data_R2','pred_perf_rsq','pred_perf_pval',
                        'nested_cv_fisher_pval', 'rho_avg', 'rho_se', 'rho_zscore', 'rho_avg_squared', 'zscore_pval',
                        'cv_rho_avg', 'cv_rho_se', 'cv_rho_avg_squared', 'cv_zscore_est', 'cv_zscore_pval', 'cv_pval_est')
write(model_summary_cols, file = model_summary_file, ncol = 26, sep = '\t')

weights_file <- output_dir %&% target_gene_id %&% '_weights_GTExv7.txt'
weights_col <- c('Geno_indicator', 'beta')
write(weights_col, file = weights_file, ncol = 2, sep = '\t')
out_beta=cbind(c(as.vector(geno_indicator)),matrix(NA,length(as.vector(geno_indicator)),1))

weights_file_train <- output_dir %&% target_gene_id %&% '_weights_train_GTExv7.txt'
weights_col_train <- c('Geno_indicator', 'beta')
write(weights_col_train, file = weights_file_train, ncol = 2, sep = '\t')
out_beta_train=cbind(c(as.vector(geno_indicator)),matrix(NA,length(as.vector(geno_indicator)),1))

perf_measures <- nested_cv_elastic_net_perf(cis_gt, adj_expression, n_samples, n_train_test_folds, n_folds, alpha)
R2_avg <- perf_measures$R2_avg
R2_sd <- perf_measures$R2_sd
pval_est <- perf_measures$pval_est
rho_avg <- perf_measures$rho_avg
rho_se <- perf_measures$rho_se
rho_zscore <- perf_measures$rho_zscore
rho_avg_squared <- perf_measures$rho_avg_squared
zscore_pval <- perf_measures$zscore_pval

# Using the weight trained on X_train to predict all data and calculate the R2
best_fit <- perf_measures$fit
all_expr_pred <- predict(best_fit, as.matrix(cis_gt), s = 'lambda.min')
all_data_R2 <- calc_R2(adj_expression, all_expr_pred)

# Fit on all data
cv_fold_ids <- generate_fold_ids(length(adj_expression), n_folds)
fit <- tryCatch(cv.glmnet(cis_gt, adj_expression, nfolds = n_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, keep = TRUE),
                error = function(cond) {message('Error'); message(geterrmessage()); list()})
if (length(fit) > 0) {
  cv_R2_folds <- rep(0, n_folds)
  cv_corr_folds <- rep(0, n_folds)
  cv_zscore_folds <- rep(0, n_folds)
  cv_pval_folds <- rep(0, n_folds)
  best_lam_ind <- which.min(fit$cvm)
  for (j in 1:n_folds) {
    fold_idxs <- which(cv_fold_ids == j)
    adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
    cv_R2_folds[j] <- calc_R2(adj_expression[fold_idxs], adj_expr_fold_pred)
    cv_corr_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor(adj_expr_fold_pred, adj_expression[fold_idxs]), 0)
    cv_zscore_folds[j] <- atanh(cv_corr_folds[j])*sqrt(length(adj_expression[fold_idxs]) - 3) # Fisher transformation
    cv_pval_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor.test(adj_expr_fold_pred, adj_expression[fold_idxs])$p.value, runif(1))
  }
  cv_R2_avg <- mean(cv_R2_folds)
  cv_R2_sd <- sd(cv_R2_folds)
  adj_expr_pred <- predict(fit, as.matrix(cis_gt), s = 'lambda.min')
  training_R2 <- calc_R2(adj_expression, adj_expr_pred)
  
  cv_rho_avg <- mean(cv_corr_folds)
  cv_rho_se <- sd(cv_corr_folds)
  cv_rho_avg_squared <- cv_rho_avg**2
  # Stouffer's method for combining z scores.
  cv_zscore_est <- sum(cv_zscore_folds) / sqrt(n_folds)
  cv_zscore_pval <- 2*pnorm(abs(cv_zscore_est), lower.tail = FALSE)
  cv_pval_est <- pchisq(-2 * sum(log(cv_pval_folds)), 2*n_folds, lower.tail = F)

  # Old way
  pred_perf <- summary(lm(adj_expression ~ fit$fit.preval[,best_lam_ind]))
  pred_perf_rsq <- pred_perf$r.squared
  pred_perf_pval <- pred_perf$coef[2,4]
  
  if (fit$nzero[best_lam_ind] > 0) {  
    out_beta[1:length(geno_indicator),2]=fit$glmnet.fit$beta[,best_lam_ind]
    write.table(out_beta,file=weights_file,sep="\t",row.names = FALSE,col.names =FALSE,quote=F )

out_beta_train[1:length(geno_indicator),2]=best_fit$glmnet.fit$beta[,perf_measures$best_lambda]
    write.table(out_beta_train,file=weights_file_train,sep="\t",row.names = FALSE,col.names =FALSE,quote=F )
        
    model_summary <- c(target_gene_id, alpha, ncol(cis_gt), fit$nzero[best_lam_ind],best_fit$nzero[perf_measures$best_lambda], fit$lambda[best_lam_ind], R2_avg, R2_sd, cv_R2_avg, cv_R2_sd, training_R2, all_data_R2, pred_perf_rsq,pred_perf_pval,pval_est,
                       rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval, cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
  } else {
    model_summary <- c(target_gene_id, alpha, ncol(cis_gt), 0, best_fit$nzero[perf_measures$best_lambda], fit$lambda[best_lam_ind], R2_avg, R2_sd,
                       cv_R2_avg, cv_R2_sd, training_R2, all_data_R2, pred_perf_rsq,pred_perf_pval,pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
                       cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
  }
} else {
  model_summary <- c(target_gene_id,  alpha, ncol(cis_gt), 0,0, NA, R2_avg, R2_sd, NA, NA, NA, NA,NA,NA, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
                     NA, NA, NA, NA, NA, NA)
}
write(model_summary, file = model_summary_file, append = TRUE, ncol = 26, sep = '\t')
