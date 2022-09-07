import sys

trait=sys.argv[1]
curr_index=sys.argv[2]

sum_file=open("/work/long_lab/jingni/project/TF_TWAS/Simulation/Simulation_Revision/Sim_GWAS/Sim_TF/summary/GTEx_Pheno"+trait+"_"+curr_index+"_model_summaries_R0.01.txt",'w')
we_file=open("/work/long_lab/jingni/project/TF_TWAS/Simulation/Simulation_Revision/Sim_GWAS/Sim_TF/summary/GTEx_Pheno"+trait+"_"+curr_index+"_weights.txt",'w')
cov_file=open("/work/long_lab/jingni/project/TF_TWAS/Simulation/Simulation_Revision/Sim_GWAS/Sim_TF/summary/GTEx_Pheno"+trait+"_"+curr_index+"_covariances.txt",'w')
for chrom in range(1,23):
    curr_chrom=str(chrom)
    sum_chr=open("/work/long_lab/jingni/project/TF_TWAS/Simulation/Simulation_Revision/Sim_GWAS/Sim_TF/summary/GTEx_Pheno"+trait+"_"+curr_index+"_chr"+curr_chrom+"_model_summaries.txt")
    head_sum=sum_chr.readline()
    if curr_chrom == "1":
        sum_file.write("gene,genename,pred.perf.R2,n.snps.in.model,pred.perf.pval,pred.perf.qval")
    sel_gene_list=[]
    for line in sum_chr:
        if (line.split("\t")[9].strip())!="NA":
            if float(line.split("\t")[9].strip())>0.01:
                sel_gene_list.append(line.split("\t")[0].strip())
                sum_file.write("\n"+line.split("\t")[0].strip()+","+line.split("\t")[1].strip()+","+line.split("\t")[9].strip()+","+line.split("\t")[6].strip()+","+line.split("\t")[11].strip()+",NA")
    we_chr=open("/work/long_lab/jingni/project/TF_TWAS/Simulation/Simulation_Revision/Sim_GWAS/Sim_TF/summary/GTEx_Pheno"+trait+"_"+curr_index+"_chr"+curr_chrom+"_weights.txt")
    head_we=we_chr.readline()
    if curr_chrom=="1":
        we_file.write("rsid,gene,weight,ref_allele,eff_allele")
    for line in we_chr:
        if line.split("\t")[0].strip() in sel_gene_list:
            we_file.write("\n"+line.split("\t")[1].strip()+","+line.split("\t")[0].strip()+","+line.split("\t")[5].strip()+","+line.split("\t")[3].strip()+","+line.split("\t")[4].strip())
    cov_chr=open("/work/long_lab/jingni/project/TF_TWAS/Simulation/Simulation_Revision/Sim_GWAS/Sim_TF/summary/GTEx_Pheno"+trait+"_"+curr_index+"_chr"+curr_chrom+"_covariances.txt")
    head_cov=cov_chr.readline()
    if curr_chrom=="1":
        cov_file.write("GENE RSID1 RSID2 VALUE")
    curr_list=[]
    for line in cov_chr:
        if line.split(" ")[0].strip() in sel_gene_list:
            pos_1=line.split(" ")[1].strip().split("_")[0].strip()+"_"+line.split(" ")[1].strip().split("_")[1].strip()
            pos_2=line.split(" ")[2].strip().split("_")[0].strip()+"_"+line.split(" ")[2].strip().split("_")[1].strip()
            rs = pos_1+" "+pos_2+" "+line.split(" ")[0].strip()
            rs2 = pos_2+" "+pos_1+" "+line.split(" ")[0].strip()
            if (rs not in curr_list) and (rs2 not in curr_list):
                curr_list.append(rs)
                curr_list.append(rs2)
                cov_file.write("\n"+line.split(" ")[0].strip()+" "+pos_1+" "+pos_2+" "+line.split(" ")[3].strip())
