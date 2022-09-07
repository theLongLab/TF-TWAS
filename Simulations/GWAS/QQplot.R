args <- commandArgs(trailingOnly=T)

# input matrix
index<-args[1]


qqPlot <- function(pval, truncate = FALSE, ylim=NULL, thinThreshold=NULL, ci=TRUE, ...) 
{
  # pvalue is a vector of p-values
  # truncate is T/F whether to truncate the y-axis (observed) to the same limits as the x-axis (expected)
  # thin is whether to thin insignificant p-values
  
  pval <- -log10(sort(pval)) # sort() removes NAs
  n <- length(pval)
  a <- 1:n
  b <- a/n
  x <- -log10(b)
  
  if (!is.null(thinThreshold)){
    breaks <- seq(from=0, to=thinThreshold, length.out=11)
    pval.cut <- cut(pval, breaks=breaks, right=FALSE)
    #quant <- quantile(pval[pval < thinThreshold], probs=1:10/10)
    #pval.cut <- cut(pval, breaks=c(0, quant), right=FALSE)
    ind.list <- list()
    for (level in levels(pval.cut)){
      sel <- which(pval.cut == level)
      ind.list[[level]] <- sel[sample.int(length(sel), min(1000, length(sel)))]
    }
    ind.list[["max"]] <- which(pval >= thinThreshold)
    ind <- unlist(ind.list, use.names=FALSE)
    ind <- sort(ind) # sorting necessary for polygon, below
  } else {
    ind <- 1:n
  }
  
  
  char <- rep(1,n)
  if(!is.logical(truncate) | truncate){
    if (is.logical(truncate)){
      maxx <- max(x)+2  
    } else {
      maxx <- min(truncate, max(pval))
    }
    
    ylm <- c(0,maxx)
    ylb <- expression(paste(-log[10], "(observed P) - truncated"))
    nx <- length(which(pval > maxx))
    if(nx > 0){
      pval[1:nx] <- maxx
      char[1:nx] <- 2
    }
    
  } else {
    ylm <- ylim
    ylb <- expression(paste(-log[10], "(observed P)"))
  }
  plot(x[ind], pval[ind], type = "n", ylim = ylm, ylab = ylb,
       xlab = expression(paste(-log[10], "(expected P)")), ...)
  # upper and lower have already been subset
  if (ci){
    upper <- qbeta(0.025, a[ind], rev(a)[ind])
    lower <- qbeta(0.975, a[ind], rev(a)[ind])
    
    polygon(-log10(c(b[ind], rev(b[ind]))), -log10(c(upper, rev(lower))), density=NA, col="gray")
  }
  points(x[ind], pval[ind], pch = char, ...)  
  abline(0,1,col="red")  
}

filename=paste("/work/long_lab/jingni/project/TF_TWAS/Simulation/Simulation_Revision/Sim_GWAS/Sim_SummStat/pheno",index,"/ALL.GRCh38.overlap_GTEx.clean.assoc.logistic",sep="")
rf<-read.table(file=filename,header = TRUE, sep=",")
df<-as.data.frame(rf)

pvals <-df$P

#outfile=paste("/work/long_lab/jingni/project/TF_TWAS/Simulation/Simulation_Revision/Sim_GWAS/Sim_SummStat/pheno",index,"/QQplot.pdf",sep="")
#pdf(file = outfile,   # The directory you want to save the file in
#    width = 10, 
#    height = 8)
#qqPlot(pvals)
#qqPlot(pvals, thinThreshold=2)
#qqPlot(pvals, truncate=TRUE)
#qqPlot(pvals, truncate=10)
#dev.off()

figure_name=paste("/work/long_lab/jingni/project/TF_TWAS/Simulation/Simulation_Revision/Sim_GWAS/Sim_SummStat/pheno",index,"/QQplot.png",sep="")
png(figure_name,width=5,height=4,units = 'in', res = 600)
qqPlot(pvals)
qqPlot(pvals, thinThreshold=2)
qqPlot(pvals, truncate=TRUE)
qqPlot(pvals, truncate=10)
dev.off()
