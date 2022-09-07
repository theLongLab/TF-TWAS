import sys
import random

index=sys.argv[1]

num_tfvar=50000 # randomly select 50K

asso_file=open("/path/to/ALL.GRCh38.overlap_GTEx_beta.assoc.logistic")
chrpos_list=[]
header=asso_file.readline()
for line in asso_file:
    curr_chrpos=line.split(",")[0].strip()+"_"+line.split(",")[2].strip()
    chrpos_list.append(curr_chrpos)
asso_file.close()

random_chrpos_list=random.sample(chrpos_list, num_tfvar)
#print(random_chrpos_list)

chrpos_rs_dict={}
chrpos_snp_dict={}
out_snp=open("/out_dir/SNP_Anno_50k_"+index+".txt",'w')
out_snp.write("SNP\tvarID\tchr\tpos\tref\teffect")
for chrom in range(1,23):
    curr_chrom=str(chrom)
    bim_file=open("/Path/to/Overlap_HG_GTEx_chr"+curr_chrom+".txt") # The overlap snp between GTEx and 1000HG
    header=bim_file.readline()
    for line in bim_file:
        curr_chrpos=line.split("\t")[0].strip()
        if curr_chrpos in random_chrpos_list:
            out_snp.write("\n"+line.strip())
        curr_snp=line.split("\t")[1].strip()
        curr_rs=line.split("\t")[0].strip()
        chrpos_rs_dict[curr_chrpos]=curr_rs
        chrpos_snp_dict[curr_chrpos]=curr_snp

    bim_file.close()

asso_file2=open("/path/to/ALL.GRCh38.overlap_GTEx_beta.assoc.logistic")
header=asso_file2.readline()
output=open("/our_dir/TF_Binding_b38_"+index+".txt",'w')
output.write("SNP_ID\tchr\tposition\tA1\tA2\tEA\tEAF\tN_Total\tBeta\tSE\tP\tTF0")
for line in asso_file2:
    curr_chrpos=line.split(",")[0].strip()+"_"+line.split(",")[2].strip()
    curr_rs=chrpos_rs_dict[curr_chrpos]
    curr_snp=chrpos_snp_dict[curr_chrpos]
    curr_ref=curr_snp.split("_")[2].strip()
    curr_alt=curr_snp.split("_")[3].strip()
    curr_chr=line.split(",")[0].strip()
    curr_pos=line.split(",")[2].strip()
    curr_rs=curr_chrpos
    curr_ea=str(0.01)
    curr_eaf="0.01"
    curr_N="2548"
    curr_beta=str(line.split(",")[6].strip())
    curr_se=line.split(",")[7].strip()
    curr_p=line.split(",")[11].strip()
    if curr_chrpos in random_chrpos_list:
        curr_tf="1"
    else:
        curr_tf="0"
    output.write("\n"+curr_rs+"\t"+curr_chr+"\t"+curr_pos+"\t"+curr_ref+"\t"+curr_alt+"\t"+curr_ea+"\t"+curr_eaf+"\t"+curr_N+"\t"+curr_beta+"\t"+curr_se+"\t"+curr_p+"\t"+curr_tf)


