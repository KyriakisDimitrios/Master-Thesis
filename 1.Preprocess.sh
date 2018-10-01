# ================================ DIRECTORIES ============================== #
WORKDIR  = "/home/kyriakis/Master_Thesis/"
GENDIR   = WORKDIR + "/0_Data/Genome"
DDRDIR   = WORKDIR + "/0_Data/Sample_Saurata-ddRAD-GR/"
FASTQDIR = WORKDIR + "/0_Data/FastQC_10_10_2017/"
RADDIR   = WORKDIR + "/1_Process_Tags/"  
SAMDIR   = WORKDIR + "/2_Align_Bowtie/"
RMAPDIR  = WORKDIR + "/3_Ref_Map/"
# =========================================================================== #


# ===================== Checking the quality of the Reads with FASTQC ====================================== #
fastqc -t 4 -o $FASTQDIR $DDRDIR/Saurata-ddRAD-GR_NoIndex_L004_R1_001.fastq.gz

# Citation
# Andrews, S. (2014). Fastqc a quality control tool for high throughput sequence data
# =========================================================================== #


# =============================== PROCESS RADTAGS =============================================#
# --------------------------------------------------------------------------------------------- #
# The first requirement is to demultiplex, or sort, the raw data to recover the individual
# samples in the Illumina library. While doing this, we will use the Phred scores provided in 
# the FASTQ files to discard sequencing reads of low quality.

# The process_radtags is a function of STACKS v.1.46 software

process_radtags 
	-P                        #  Files contained within directory specified by '-p' are paired.
	-i gzfastq                #  Input file type

	-1 ./data/R1.fastq.gz     #  First input file in a set of paired-end sequences.
	-2 ./data/R2.fastq.gz     #
	-b barcode                #  Path to a file containing barcodes for this run.
	
	--inline_inline           #  Barcode occurs in FASTQ header
	--renz_1 sbfI             #  The restriction enzyme used 
	--renz_2 sphI             #  The second restriction enzyme used (cut site occurs on he 		
							  #  paired-end read)
	
	-c                        #  Clean data, remove any read with an uncalled base.
	-q                        #  Discard reads with low quality scores.
	-r                        #  Rescue barcodes and RAD-Tags.
	-t 100                    #  Truncate final read length to this value.
	-o ./new_version          #  path to output the processed files.
	-D
#
# Others:
# p — path to a directory of files.
# E — specify how quality scores are encoded, 'phred33' (Illumina 1.8+, Sanger, default) or 
#             'phred64' (Illumina 1.3 - 1.5).       ******************
# D — capture discarded reads to a file.

# Citation
# Catchen, J. M. (2013). Stacks: an analysis tool set for population genomics. Molecular ecology 22,
# 3124–3140. doi:10.1111/mec.12354.Stacks
# --------------------------------------------------------------------------------------------- #

