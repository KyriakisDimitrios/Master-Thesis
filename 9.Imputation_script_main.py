# -*- coding: utf-8 -*-
"""
Created on 2017
Project: Master Thesis
Scanning of genetic variants and genetic mapping of phenotypic traits in gilthead seabream through
ddRAD sequencing data analysis

@author: Dimitris Kyriakis
"""

#========================== LOAD Library  ==========================#
from Imputation_Functions import *
#-------------------------------------------------------------------#

#=== ENABLE MODES ====#
roc         = 1 # PLot roc curves
#--------------------#

#================ Choose Variant and  Classifier ================== #
# clas = "1"
# name = "KNN"
# clas = "2"
# name = "LG"
clas="4"
name = "SVM"
# clas="3"
# name = "RF"
#------------------------------------------------------------------ #


#================ Choose Variants for analysis ==================== #
neigh   = 0 # KEEP ONLY NEIGHBOUR SNPS +-40
ld  = 1 # If 1 => Only variants in same chrom 
rev = 0	  # If 1 => Only variants in different chrom

#-------------------------------------------------------------------#
# print("\n\n"+"Classifier = "+name)
# print("Only the same Chromosome = {}".format(ld))
output2 = open("Impu_Results_LD{}.txt".format(ld),"w")
#-------------------------------------------------------------------#

#=============================================================================#
#============================ CHOOSE DATA SET ================================#
#=============================================================================#
Train = "../FullDataset_105.txt"
Test = "../Impute_Dataset.txt"
HomR = 0
Hetr = 1
HomA = 2
#-------------------------------------------------------------------#

# ==================== READ TRAIN DATA ============================ #
# With no missing
Data = np.genfromtxt(Train, delimiter='\t', dtype=str).T 
# Chrom name of Variants
Var_Names = np.asarray([x[:x.index(":")] for x in Data[0,1:].tolist()] ) #
# All non missing Variants : Rows = Samples, Col = Variants
Full_Data = (Data[1:,1:]).astype(np.float)
# Names of the samples
Samples_Name = Data[1:,0].reshape(Data.shape[0]-1,1)
# Names of the Features
Features = Data[0,1:]
Features_names = [i[:-4] for i in Features]
#--------------------------------------------------------------------#
counter =0
Variant_num = 11
header = Features.flatten().tolist()


for Variant_num in range(1,3028):

	print(Variant_num)
	# ====================== READ TEST DATA ============================ #
	# Name of Variant to Test (with missing value(s))
	SNP = np.genfromtxt(Test, delimiter='\t',dtype=str).T[0,Variant_num]
	# Genotype of the SNP
	Test_Data_raw = np.genfromtxt(Test, delimiter='\t',missing_values='-').T[1:,Variant_num]
	# Index of missing values
	missing_index = np.argwhere(np.isnan(Test_Data_raw)).flatten().tolist()
	# Name of missing Samples
	Test_Data_samples = Samples_Name[missing_index].flatten().tolist()
	#--------------------------------------------------------------------#


	# ======================= Prepare Data sets ======================== #
	# Index of Variants in the same Chrom
	Keep_Vars  = np.where(Var_Names==SNP[:SNP.index(":")])[0].tolist()
	# Delete the samples with missing values for the SNP in Train Data
	Train_Data = np.delete(Full_Data,missing_index,0)
	# The Samples Data of the missing 
	New_Data = Full_Data[missing_index,:].reshape(len(missing_index),-1)
	
	
	#============ IF rev=1 => Keep only Vars in different chrom
	if rev==1:
		Train_Data = np.delete(Train_Data,Keep_Vars,1)
		New_Data = np.delete(New_Data,Keep_Vars,1)
	
	# ============ KEEP ONLY NEIGHBOUR SNPS +-40 ==============#
	if neigh ==1:
		pos_list = temp_list = [ int(x.split(":")[1]) for x in Features[Keep_Vars]]
		snp_pos = int(SNP.split(":")[1])
		temp_list.append(snp_pos) ; temp_list.sort() ; pos_index = temp_list.index(snp_pos)
		Train_Data = Train_Data[:,Keep_Vars]
		num_features = Train_Data.shape[1]
		if pos_index<40:
			start = 0 
		else:
			start=pos_index-40
		if num_features - pos_index>= 0:
			end = 40 + pos_index
		else:
			end = num_features
		Train_Data = Train_Data[:,start:end]
		New_Data = New_Data[:,Keep_Vars]
		New_Data = New_Data[:,start:end]
	#-----------------------------------------------------------#


    #============ IF ld=1 => Keep only Vars in same chrom
	if ld==1:
		if len(Keep_Vars):
			Train_Data = Train_Data[:,Keep_Vars]
			New_Data = New_Data[:,Keep_Vars]




	# Known Labels for theb SNP that will Impute
	Raw_Labels     =  np.delete(Test_Data_raw,missing_index)
	Known_Labels   = Raw_Labels.astype(int).astype(str)
	Known_Labels_ROC   = Raw_Labels.astype(int)
	#--------------------------------------------------------------------#


	#======================= COUNT GENOTYPES =============================#
	Num_HomR = sum(Raw_Labels == 0.0)
	Num_Hetr = sum(Raw_Labels == 1.0)
	Num_HomA = sum(Raw_Labels == 2.0)

	if (Num_HomR==Num_Hetr==0) or (Num_HomR==Num_HomA==0) or (Num_HomA==Num_Hetr==0):
			continue

	# print("Number of Homogygous Reference : {}".format(Num_HomR))
	# print("Number of Heterogygous : {}".format(Num_Hetr))
	# print("Number of Homogygous  Alternative: {}".format(Num_HomA))
	# print("\n")

	# Set minimum Folds Based on classes
	Min_Folds = min(Num_HomR,Num_Hetr,Num_HomA)
	if Num_HomA == 0:
	    Min_Folds = min(Num_HomR,Num_Hetr)
	if Num_HomR == 0:
	    Min_Folds = min(Num_Hetr,Num_HomA)
	if Num_Hetr == 0:
	    Min_Folds = min(Num_HomR,Num_HomA)
	#---------------------------------------------------------------------#
	if Min_Folds ==2 or Min_Folds ==1:
		continue
	if Min_Folds >=5:
		Min_Folds = 5
	
	#=================== CONVERT ==================================\
	X = Train_Data.copy()
	New = New_Data.copy()

	
	#========================== Cross Validation =========================#
	mean_acc_f,pipe,predict = CV_Kyriakis(SNP,X,Known_Labels,New,name,Min_Folds,Test_Data_samples)
	# ============================= ROC - AUC =========================== # 
	if roc ==1 and mean_acc_f >= 0.85:
		uni_clas  =np.unique(Known_Labels_ROC) 
		if len(uni_clas) ==2:
			micro=0	
			if 1 in uni_clas:
				mean_auc_0 = ROC_Kyriakis(SNP,X,Known_Labels_ROC,name,mean_acc_f,predict,ld,Min_Folds,pipe,1)
			else:
				mean_auc_0 = ROC_Kyriakis(SNP,X,Known_Labels_ROC,name,mean_acc_f,predict,ld,Min_Folds,pipe,2)	
		else:
			micro,macro = ROC_VS(SNP,Train_Data,Known_Labels_ROC,name,mean_acc_f,predict,ld,Min_Folds,pipe)
		if (micro>=0.9 and macro>=0.9) or mean_auc_0>=0.9:
			for i in range(len(predict['Id'])):
				a = np.where(Samples_Name ==predict['Id'][i])
				Test_Data_raw[a[0]] = predict['Predicted'][i]
			if  Variant_num==1:
				Imputed_Matrix =np.concatenate((Full_Data,Test_Data_raw.reshape(105,1)),axis =1)
			else:
				Imputed_Matrix =np.concatenate((Imputed_Matrix,Test_Data_raw.reshape(105,1)),axis =1)
			output2.write(SNP+"\t"+str(predict)+"\n")
			header.append(SNP)


import pandas as pd 
df = pd.DataFrame(Imputed_Matrix.astype(int))
df.to_csv("Data_With_Imputation.tsv",index=False,header=header,sep='\t')
output2.close()
#---------------------------------------------------------------------#