# ======================== Align data against a reference genome ============================= #
# --------------------------------------------------------------------------------------------- #
# Align these samples using a standard alignment program,  Bowtie2.
# Creating an index/database of your reference (this speeds up the mapping).

# ------------------------------------ BOWTIE2 -------------------------------------------------- #
# Bowtie2 v.2.3.0

bowtie = "bowtie2 -S -p 12 --very-sensitive -x {0}".format(GenomeDIR+'Bowtie_Build/sparus_aurata')

# Align reads
os.system('bowtie2-build {0}/Saurata_genome.fa {0}/Bowtie_Build/sparus_aurata'.format(GENDIR))

bowtie_params ="bowtie2 -p 20 --end-to-end --sensitive --no-unal"

for name in Barcodes:
	# ==== Merge PARENTS ==== #
	if "Br" in name: 
		os.system("zcat {1}/{0}_?.1.fq.gz|gzip  > {1}/{0}.1.fq.gz".format(name,RADDIR))
		os.system("zcat {1}/{0}_?.2.fq.gz|gzip  > {1}/{0}.2.fq.gz".format(name,RADDIR))
		os.system("zcat {1}/{0}_?.rem.1.fq.gz|gzip  > {1}/{0}.rem.1.fq.gz".format(name,RADDIR))
		os.system("zcat {1}/{0}_?.rem.2.fq.gz|gzip  > {1}/{0}.rem.2.fq.gz".format(name,RADDIR))
	os.system("echo {0} >> {1}/0_Scripts/2_SAM.out".format(name,WORKDIR))

	# =============================== ALIGN ================================= #
	os.system("{4} --rg-id {2} -x {0}/Bowtie_Build/sparus_aurata -1 {1}/{2}.1.fq.gz -2 {1}/{2}.2.fq.gz -U {1}/{2}.rem.1.fq.gz,{1}/{2}.rem.2.fq.gz -S {3}/{2}.sam".format(GENDIR,RADDIR,name,SAMDIR,bowtie_params))

	# ============================= FILTERING =============================== #
	# Save header
	os.system("samtools view -H {1}/{0}.sam > {1}/header.sam".format(name,SAMDIR))
	# Remove multi aligned | Remove >mismatches| Create Bam file
	os.system("grep -v 'XS:'|wc -l {1}/{0}.sam |paste {1}/mul_head - >> {2}/0_Scripts/2_SAM.out".format(name,SAMDIR,WORKDIR))
	os.system("grep 'XM:i:[0-3]'|wc -l {1}/{0}.sam | paste {1}/mis_head - >> {2}/0_Scripts/2_SAM.out".format(name,SAMDIR,WORKDIR))
	os.system("samtools view -F 4 -q 20 {1}/{0}.sam | grep -v 'XS:'| grep 'XM:i:[0-3]' | cat {1}/header.sam - |samtools view -b -o {1}/{0}.bam".format(name,SAMDIR))
	# Sort Bam file
	os.system("samtools sort {1}/{0}.bam -o {1}/{0}.sortred.bam".format(name,SAMDIR))
	# Delete header
	os.system("rm {1}/header.sam {1}/{0}.bam {1}/{0}.sam ".format(name,SAMDIR))# --end-to-end       entire read must align; no clipping (on)

# --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
#  -D <int>           give up extending after <int> failed extends in a row (15)
#  -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)
#  -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)
#  -L <int>           length of seed substrings; must be >3, <32 (22)
# 
# -v  2  # no more than 2 missmatches
# -n -l # max number of mismatches in the high quality “seed”, which is the the first l
# base pairs of a read

# Citation
# Langmead, B. and Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods 9,
# 357–359. doi:10.1038/nmeth.1923
# --------------------------------------------------------------------------------------------- #
