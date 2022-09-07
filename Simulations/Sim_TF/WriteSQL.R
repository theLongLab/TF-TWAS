library(sqldf)

"%&%" <- function(a,b) paste(a,b, sep = "")

args <- commandArgs(trailingOnly=T)
  
# input matrix
trait<-args[1]
curr_index<-args[2]

working_dir="./summary/"

db <- dbConnect(SQLite(),dbname=working_dir%&%"GTEx_"%&%curr_index%&%".db")
dbWriteTable(conn=db,name="extra",value=working_dir%&%"GTEx_"%&%curr_index%&%"_model_summaries_R0.01.txt",row.names=FALSE,header=TRUE)
dbWriteTable(conn=db,name="weights",value=working_dir%&%"GTEx_"%&%curr_index%&%"_weights.txt",row.names=FALSE,header=TRUE)
