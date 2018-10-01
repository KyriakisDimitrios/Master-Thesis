# ================================ DIRECTORIES ============================== #
WORKDIR  = "/home/kyriakis/Master_Thesis/"
GENDIR   = WORKDIR + "/0_Data/Genome"
DDRDIR   = WORKDIR + "/0_Data/Sample_Saurata-ddRAD-GR/"
FASTQDIR = WORKDIR + "/0_Data/FastQC_10_10_2017/"
RADDIR   = WORKDIR + "/1_Process_Tags/"
SAMDIR   = WORKDIR + "/2_Align_Bowtie/"
RMAPDIR  = WORKDIR + "/3_Ref_Map/"

import os



# =============================== 1. PSTACKS ==================================== #
# The pstacks program will extract stacks that have been aligned to a reference genome
# by an aligner (Bowtie2). Pstacks compares "stacks of reads" and forms
# putative sets of loci. These sets are used in order to detect SNPs at each locus using a
# maximum likelihood framework.

command = "stacks pstacks -p 12 -o {0}/3_Pstacks -m 3 ".format(WORKDIR)
lista = []
counter = 0
for file in os.listdir(SAMDIR):
	if file.endswith("bam"):
		file = file [:file.index(".",file.index(".")+1)]
		if file not in lista:
			counter+=1
			lista.append(file)
		command += " -f {0}/{1}.bam -i {2}".format(SAMDIR ,file,counter)
os.system(command)
# =============================================================================== #



# =================================== 2. CSTACKS ==================================== #
# A SNP catalogue was built from the parents of the cross. Cstacks created a set of all
# possible alleles expected in the progeny of the cross.

command = "stacks cstacks -b 1 -p 12 -o {0}/Gen_Cstacks/ --aligned ".format(RMAPDIR)
lista = []
for file in os.listdir(RMAPDIR+"/Pstacks/Parents/"):
	file = file [:file.index(".",file.index(".")+1)]
	if file not in lista:
		lista.append(file)
		command += " -s {1}/Pstacks/Parents/{0} ".format(file,RMAPDIR)
#print(command)
os.system(command)
# -b Batch id
# -p Multithreads





# =================================== 3. SSTACKS =================================== #
# The sets of stacks that constructed by the pstacks program searched against the catalog
# produced by cstacks. All samples in the population matched against the catalog with
# sstacks.

sstacks_com = "stacks sstacks -g -p 12 -b 1 -c {0}/Gen_Cstacks/batch_1 -o {0}/Gen_Sstacks/".format(RMAPDIR) 

# KIDS
lista = []
for file in os.listdir(RMAPDIR+"/Pstacks/"):
	if file.endswith("gz"):
		file = file [:file.index(".",file.index(".")+1)]
		if file not in lista:
	        	lista.append(file)
			file_com = sstacks_com + " -s {1}/Pstacks/{0}".format(file,RMAPDIR)
			os.system(file_com)
# PARENTS
lista = []
for file in os.listdir(RMAPDIR+"/Pstacks/Parents"):
	if file.endswith("gz"):
		file = file [:file.index(".",file.index(".")+1)]
		if file not in lista:
        		lista.append(file)
			file_com = sstacks_com + " -s {1}/Pstacks/Parents/{0}".format(file,RMAPDIR)
			os.system(file_com)



# =================================== Population =================================== #
# The populations program uses the population map to determine which groupings to
# use for calculating summary statistics, such as heterozygosity.

Main = "stacks populations "
Files = " -P {0} -O {0}/Results -M {0}/popmap1 -b 1 -k ".format(POPStaDIR)
Params = " -f p_value -t 12 --structure --vcf --vcf_haplotypes --plink "

Command = Main+Files+Params
os.system(Command)
# =============================================================================== #