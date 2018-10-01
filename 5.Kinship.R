
#  Check family relationship and indicate possible pedigree errors
# 1. Principal Components analysis (PCA)
# 2. Hierarchical clustering
# 3. KING Kinship


# === Libraries === #
library(ggplot2)
library(webshot)
library(cluster)
library(pca3d)
# ================= #


# ===================== LOAD DATA ========================= #
# TSV file with rows:samples and columns: SNPs.
# 0 = Homozygous to Reference
# 1 = Heterozygous 
# 2 = Homozygous to Alternative

# Prepare Tables
mydata = read.csv("Dataset.tsv",header = TRUE,sep="\t")
mydata[ mydata == "-" ] <- NA
attach(mydata)
Pheno    = as.data.frame(mydata[3:7])
All_Features = as.data.frame(mydata[8:dim(mydata)[2]])
Names = colnames(Pheno)
IDS = paste("Family",Fam_ID)
new = cbind(IDS,All_Features)
row.names(new) <-mydata[1]$CHROM.POS.ID.REF.ALT



# ======================== PCA =========================
# == PLOT PC Variance
pca <- prcomp(All_Features)
gr <- factor(new[,1])
library(scatterplot3d)
levels(Fam_ID)
colors1 =palette(rainbow(6))
colors1[7]="coral4"
colors <- colors1[as.numeric(Fam_ID)]


pca3d(pca, group=IDS, show.ellipses=TRUE,
      ellipse.ci=0.75, show.plane=FALSE, legend="right",labels.col=colors1)
snapshotPCA3d(file="PCA3D.png")



# ====================== Hierarchical clustering =========================

library(ape)
new2=All_Features
rownames(new2) <-mydata[1]$CHROM.POS.ID.REF.ALT

di <- dist(new2, method="euclidean")
hc <- hclust(di, method="ward.D2")

pdf("Hierar.pdf")
plot(as.phylo(hc), type = "fan", tip.color = colors,
     label.offset = 1,cex=0.9)
par(mar=c(0, 0, 0, 0))
legend(-50,-80,'bottom',inset=.05, legend=sort(unique(IDS)),
         cex=0.9,  fill=colors1, xpd = NA, bty = "n",ncol=4)
dev.off()
# ========================================================================================
# ========================================================================================




# ============================== PLINK ================================== #
# PLINK v1.90p 64-bit (14 Nov 2017)              www.cog-genomics.org/plink/1.9/                                                                                                                                     
# (C) 2005-2017 Shaun Purcell, Christopher Chang   GNU General Public License v3  
# --allow -extra -chr : allow scaffolds

plink --vcf batch_1.vcf --allow -extra -chr --make-bed --out plink
plink --vcf batch_1.vcf --allow -extra -chr --recode oxford

# ======================== KING ========================== #
# (C) 2005-2017 Shaun Purcell , Christopher Chang GNU General Public License v3
# KING 2.1 - (c) 2010-2018 Wei-Min Chen
king -b plink.bed --cluster --degree 10
king -b plink.bed --kinship --degree 10

