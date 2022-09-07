library(sqldf)

"%&%" <- function(a,b) paste(a,b, sep = "")

args <- commandArgs(trailingOnly=T)
  
# input matrix
trait<-args[1]
curr_index<-args[2]

working_dir="/work/long_lab/jingni/project/TF_TWAS/Simulation/Simulation_Revision/Sim_GWAS/Sim_TF/summary/"

db <- dbConnect(SQLite(),dbname=working_dir%&%"GTEx_Pheno"%&%trait%&%"_"%&%curr_index%&%".db")
dbWriteTable(conn=db,name="extra",value=working_dir%&%"GTEx_Pheno"%&%trait%&%"_"%&%curr_index%&%"_model_summaries_R0.01.txt",row.names=FALSE,header=TRUE)
dbWriteTable(conn=db,name="weights",value=working_dir%&%"GTEx_Pheno"%&%trait%&%"_"%&%curr_index%&%"_weights.txt",row.names=FALSE,header=TRUE)
