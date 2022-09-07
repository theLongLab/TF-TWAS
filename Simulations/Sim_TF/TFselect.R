library(data.table)
library(lme4)
library(parallel)
library(fastmatch)


args <- commandArgs(trailingOnly=T)

# input matrix
index<-args[1]

mat <- paste("/path/to/out_dir/TF_Binding_b38_",index,".txt",sep="")
work_dir<-paste('/path/to/work_dir',sep="")
setwd(work_dir)
trait <- fread(mat,header=T)
names(trait)

dim(trait)
print("trait")
print(names(trait))
print(dim(trait))

# calculate t value 
trait <- within (trait, {
freq  <-as.numeric(as.character(EAF))
beta  <-as.numeric(as.character(Beta))
se <-as.numeric(as.character(SE))
tv<-beta/se
})

trait<-trait[!duplicated(trait[,2:3]),]
names(trait)
dim(trait) 


trait<-trait[!is.na(trait$tv) & trait$chr>=1 & trait$chr<=22,]
dim(trait) 

trait0<-trait[trait$tv==0,]
n0<-dim(trait0)[1]
rm(trait0)
set.seed(123)
trait$tv[trait$tv==0]<-runif(n0,-0.001,0.001)

dim(trait) 
range(trait$tv) 

trait<-trait[order(trait$chr,trait$position),]

trait$ID<-paste(trait$chr,'-',trait$position,sep='')  
trait$pc<-cut(trait$P,breaks=c(0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,
                   0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),label=1:15)
table(trait$pc)

pf<-table(trait$pc)
# estimate the expected number of SNPs in each p-values category
nx<-ceiling(10*min(pf[7:15])*c(1e-6,1e-5-1e-6,1e-4-1e-5,1e-3-1e-4,1e-2-1e-3,1e-1-1e-2,rep(0.1,9)))
idv<-vector('list',15)

set.seed(123)
for (i in 1:15) {
bdx<-trait[trait$pc==i,]
idv[[i]]<-sample(bdx$ID,nx[i],replace = TRUE)
}
table(duplicated(unlist(idv)))
trait$chosen<-ifelse(trait$ID %in% unlist(idv),1,0)
rm(idv)
print(table(trait$chosen)) #1: chosen for deflated genome, 0: enriched with highly significant SNPs

print(names(trait)[12:12])

### LD block 100kb
names(trait)[2:3] <-c("chr","position")

KB<-100000

trait$loci <- paste0(trait$chr,'_',floor(trait$position/KB))

length(unique(trait$loci)) 

rm(bdx,i,KB,n0,nx) 


#############################################################################################################################################

trait$Pcc<-ifelse(trait$P<5e-8,1,0)
dim(trait)  
trait<-trait[!is.na(trait$tv),]
table(is.na(trait$tv)) 
dim(trait)  
table(trait$Pcc)
#      0       1
#6491120    3458
table(trait$chosen)
#     0       1
# 399535 6095043

trait <- as.data.frame(trait)
print(names(trait)[12:12])
####  Association of single TF with Chi-square

print("Final_trait")
print(names(trait))
print(dim(trait))
trait1<-setDT(trait)

outr1 <- mclapply(trait1[,c(12:12)],function(x) {summary(lmer(I(trait1$tv^2)~x+(1|trait1$loci),control = lmerControl(calc.derivs = FALSE)))$coef[2,]}, mc.cores = 6)


outr11 <- data.frame(matrix(unlist(outr1), ncol = 3, byrow = TRUE))
row.names(outr11) <- names(trait)[12:12]
colnames(outr11) <- c("beta","se","tv")
out_name=paste("TFoutr_",index,".txt",sep="")
write.table(outr11,out_name,row.names=T,sep="\t",quote=F)


rf<-read.table(file=out_name,header = TRUE, sep="\t")
df<-as.data.frame(rf)
rownames=rownames(df)
print(rownames)
beta<-df$beta
se<-df$se
z<-beta/se

pval<-2*pnorm(q=z, lower.tail=FALSE)
#print(pval)

df$p<-pval

if(pval<0.05){
	#print(df)
	out_name=paste("TFoutr_P_",index,".csv",sep="")
	write.csv(df,out_name, quote=F,row.names = rownames)
}
