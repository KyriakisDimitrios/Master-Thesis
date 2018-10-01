

# Convert VCF file to tsv file. Columns: SNPs, Rows: Samples. 
# 0 = Homozygous to Reference
# 1 = Heterozygous 
# 2 = Homozygous to Alternative

# ==  Libraries == # 
import numpy as np
import re
import pandas as pd
# ================ #

Samples_ord = []
Dict = {}
Order_SNPS=[]

# == Thresholds == #
AF = 0.05
NS = 100
# == COUNTERS == #
counter =0
counter2 =0
counter3 =0

transitions = ["A_G","G_A","C_T","T_C"]
Ts_Tv = {"Ts":0,"Tv":0}

def snp_conv(x):
	geno = x[:3]
	if geno[0] == "0":
		if geno[2] == "0":
			return "0"
		else:
			return "1"	
	if geno[0] == ".":
		return "-"
	return "2"


# ======================== SNPS ========================= #
vcf = open("batch_1.vcf","r")

for line in vcf:
	# ==================== SKIP ===================== #
	if line.startswith("##"):
		continue
	line2 = line
	line = line.rstrip().split("\t")
	# =================== HEAD ====================== #
	if line[0] == "#CHROM":
		Samples_ord = np.array(line[9:])
		rm_index =  ['Br_' not in x for x in  Samples_ord ]
		Samples_ord = [x[:-8] for x in Samples_ord[rm_index]]
		arg_sort = np.argsort(Samples_ord)
		Samples_ord = np.asarray(Samples_ord)[arg_sort]
		out.write("CHROM:POS:ID:REF:ALT\t{}\n".format("\t".join(Samples_ord))) 
		continue
	# ================= MAF < Threshold ====================== #
	if float(line[7].split(";")[1].split("=")[1]) <= AF:
		counter2 +=1
		continue


	# ================ SAMPLES < THreshold =================== #
	if int(line[7].split(";")[0].split("=")[1]) < NS or int(line[7].split(";")[0].split("=")[1]) == 105:
		counter3 +=1
		continue
	# ------------------------------------------------ #

	counter +=1
	col_name = (":").join(line[:5])
	if ("_").join(line[3:5]) in transitions:
		Ts_Tv["Ts"] +=1
	else:
		Ts_Tv["Tv"] +=1
	snps = np.array(list(map(snp_conv,line[9:])))[rm_index]
	snps = snps[arg_sort]
	Order_SNPS.append(col_name)
	Dict[col_name] = snps
	out.write("{}\t{}\n".format(col_name,"\t".join(snps)))


# ==== Summary Statistics === # 
print("Total:",counter3+counter2+counter,"\n","="*10)
print("Active :   ",counter)
print("AF < {} :".format(AF),counter2)
print("NS < {} : ".format(NS),counter3)
print("Ts/Tv  =    {:.4f}".format(Ts_Tv["Ts"]/Ts_Tv["Tv"]))


# ======= CONVERT VCF TO TSV ==== #
Dict["CHROM:POS:ID:REF:ALT"]=Samples_ord
df = pd.DataFrame(data=Dict)
Order_SNPS_all = ["CHROM:POS:ID:REF:ALT"]+Phenodata+Order_SNPS
# Phenodata: Matrix with columns: phenortypes, rows: samples
df = df[Order_SNPS_all]
df = df.set_index("CHROM:POS:ID:REF:ALT")
df.to_csv("Dataset.tsv", sep='\t')
