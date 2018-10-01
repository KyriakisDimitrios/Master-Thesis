library(lme4)
library(ggplot2)
library(webshot)
library(GWASTools)


# ===================== LOAD DATA ========================= #
mydata = read.csv("Dataset.tsv",header = TRUE,sep="\t")
mydata[ mydata == "-" ] <- NA
attach(mydata)
Pheno    = as.data.frame(mydata[c(3,6,7,4)])#3:7
All_Features = as.data.frame(mydata[8:dim(mydata)[2]])
Names = colnames(Pheno)
features = colnames(All_Features)
IDS = paste("Family",Fam_ID)
# --------------------------------------------------------- #

# ====================== MIXED MODELS ==================================
# Initialize p-values, b coefficients
pvals_Fat = rep(0,length(All_Features))
pvals_Weight = rep(0,length(All_Features))
pvals_Tag_Weight = rep(0,length(All_Features))
pvals_Length_Width = rep(0,length(All_Features))

coefs_Fat = rep(0,length(All_Features))
coefs_Weight = rep(0,length(All_Features))
coefs_Tag_Weight = rep(0,length(All_Features))
coefs_Length_Width = rep(0,length(All_Features))
#---------------------------------------------

# Iterate between the SNPs
for (i in 1:length(All_Features)){
  #=== FAT
  x = All_Features[,i]
  # Constract a linear mixed Model with SNP i predicting a phenotypic trait
  # IDS = Family id 
  m1 = lmer(Fat~x+(1|IDS),REML=FALSE)
  m2 <- update(m1, .~. - x)
  pvals_Fat[i]=(anova(m1,m2)$`Pr(>Chisq)`[2])
  if (pvals_Fat[i] == 1){
  		coefs_Fat[i] = 0
  } else{
  	coefs_Fat[i] = coef(m1)[[1]]$x[1]
  }
  #=== Weight
  m3 = lmer(Weight~x+(1|IDS),REML=FALSE)
  m4 <- update(m3, .~. - x)
  pvals_Weight[i]=(anova(m3,m4)$`Pr(>Chisq)`[2])
  if (pvals_Weight[i] == 1){
  		coefs_Weight[i] = 0
  } else{
  	coefs_Weight[i] = coef(m3)[[1]]$x[1]
  }
  #=== Tag Weight
  m5 = lmer(Tag_weight~x+(1|IDS),REML=FALSE)
  m6 <- update(m5, .~. - x)
  pvals_Tag_Weight[i]=(anova(m5,m6)$`Pr(>Chisq)`[2])
  if (pvals_Tag_Weight[i] == 1){
  		coefs_Tag_Weight[i] = 0
  } else{
  	coefs_Tag_Weight[i] = coef(m5)[[1]]$x[1]
  }
  #=== Length_Width
  m7 = lmer(Length_Width~x+(1|IDS),REML=FALSE)
  m8 <- update(m7, .~. - x)
  pvals_Length_Width[i]=(anova(m7,m8)$`Pr(>Chisq)`[2])
  if (pvals_Length_Width[i] == 1){
  		coefs_Length_Width[i] = 0
  } else{
  	coefs_Length_Width[i] = coef(m7)[[1]]$x[1]
  }
}



# ===== Chrom names
temp1=strsplit(colnames(All_Features),"\\.")
chroms = as.data.frame(matrix(unlist(temp1), ncol=5, byrow=TRUE))$V1
chr1 = which(grepl("chr1$",chroms))
chr2 = which(grepl("chr2$",chroms))
chr3 = which(grepl("chr3$",chroms))
chr4 = which(grepl("chr4$",chroms))
chr5 = which(grepl("chr5$",chroms))
chr6 = which(grepl("chr6$",chroms))
chr7 = which(grepl("chr7$",chroms))
chr8 = which(grepl("chr8$",chroms))
chr9 = which(grepl("chr9$",chroms))
chr10 = which(grepl("chr10",chroms))
chr11 = which(grepl("chr11",chroms))
chr12 = which(grepl("chr12",chroms))
chr13 = which(grepl("chr13",chroms))
chr14 = which(grepl("chr14",chroms))
chr15 = which(grepl("chr15",chroms))
chr16 = which(grepl("chr16",chroms))
chr17 = which(grepl("chr17",chroms))
chr18 = which(grepl("chr18",chroms))
chr19 = which(grepl("chr19",chroms))
chr20 = which(grepl("chr20",chroms))
chr21 = which(grepl("chr21",chroms))
chr22 = which(grepl("chr22",chroms))
chr23 = which(grepl("chr23",chroms))
chr24 = which(grepl("chr24",chroms))

jindx = c(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24)
jindx = chroms[c(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24)]

jindx = grepl("chr",chroms)
dindx = grepl("scaf",chroms)

# length(jindx)
chroms
ordered_index = order(as.integer(substring(chroms[jindx], first = 4)))



# Independent SNPs based on LD using plink = 497
Threshold = 0.05/497

#=== FAT
png(filename="Fat_QQplot.png", units="px", width=1600, height=1600, res=300)
qqPlot(pvals_Fat)
dev.off()

png(filename="Fat_Manhattan.png",  units="px", width=1600, height=1600, res=300)
manhattanPlot(pvals_Fat[jindx][ordered_index],chroms[jindx][ordered_index],signif=Threshold,ylim = c(0,5))
dev.off()

png(filename="Fat_Manhattan_Scaf.png",  units="px", width=1600, height=1600, res=300)
manhattanPlot(pvals_Fat[dindx],chroms[dindx],signif=Threshold,ylim = c(0,5))
dev.off()


# #=== Weight
png(filename="Weight_QQplot.png",  units="px", width=1600, height=1600, res=300)
qqPlot(pvals_Weight)
dev.off()

png(filename="Weight_Manhattan.png",  units="px", width=1600, height=1600, res=300)
manhattanPlot(pvals_Weight[jindx][ordered_index],chroms[jindx][ordered_index],signif=Threshold,ylim = c(0,5))
dev.off()

png(filename="Weight_Manhattan_Scaf.png",  units="px", width=1600, height=1600, res=300)
manhattanPlot(pvals_Weight[dindx],chroms[dindx],signif=Threshold,ylim = c(0,5))
dev.off()






#=== Tag_Weight
png(filename="Tag_QQplot.png",  units="px", width=1600, height=1600, res=300)
qqPlot(pvals_Tag_Weight)
dev.off()

png(filename="Tag_Manhattan.png",  units="px", width=1600, height=1600, res=300)
manhattanPlot(pvals_Tag_Weight[jindx][ordered_index],chroms[jindx][ordered_index],signif=Threshold,ylim = c(0,5))
dev.off()

png(filename="Tag_Manhattan_Scaf.png",  units="px", width=1600, height=1600, res=300)
manhattanPlot(pvals_Tag_Weight[dindx],chroms[dindx],signif=Threshold,ylim = c(0,5))
dev.off()





#=== Length_Width
png(filename="Len_Wid_QQplot.png",  units="px", width=1600, height=1600, res=300)
qqPlot(pvals_Length_Width)
dev.off()

png(filename="Len_Wid_Manhattan.png",  units="px", width=1600, height=1600, res=300)
manhattanPlot(pvals_Length_Width[jindx][ordered_index],chroms[jindx][ordered_index],signif=Threshold,ylim = c(0,5))
dev.off()

png(filename="Len_Wid_Manhattan_Scaf.png",  units="px", width=1600, height=1600, res=300)
manhattanPlot(pvals_Length_Width[dindx],chroms[dindx],signif=Threshold,ylim = c(0,5))
dev.off()



