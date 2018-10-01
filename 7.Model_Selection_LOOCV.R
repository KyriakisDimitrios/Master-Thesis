library(lme4)
library(MXM)


source("./.../Cross_Validation_LOOVC_Functions.R")

# ----------------------------------------------------------------------------- #

# =================== LOAD GENOTYPES WITH 10% MISSING ========================= #
mydata1 = read.csv("Dataset.tsv",header = TRUE,sep="\t")
Most_Freq_Impute="y"
# If there are missing values in TSV turn on

# Initializations
z <- mydata
Fat    <- z[, 3]
Length_Width <- z[, 4]
Tag_Weight <- z[, 6]
Weight <- z[, 7]


id <- z[, 2]
x <- as.matrix(z[, -c(1:7)])
N <- length(Fat)
er1 <- er0 <- numeric(N)
m <- length( unique(id) )


# ======================== INITIAL DATAFRAME WITH RESULTS ====================== #
Pheno_Names=c("Fat","Weight","Tag_Weight","Length_Width")
results_df <- data.frame(matrix(ncol = 4, nrow = 22))
algorithms <- c("Mean", "Median",
	"Mean_SES_001_2","Median_SES_001_2",
	"Mean_SES_001_3","Median_SES_001_3",
	"Mean_SES_001_4","Median_SES_001_4",
	"Mean_SES_001_5","Median_SES_001_5",
	"Mean_SES_005_2","Median_SES_005_2",
	"Mean_SES_005_3","Median_SES_005_3",
	"Mean_SES_005_4","Median_SES_005_4",
	"Mean_SES_005_5","Median_SES_005_5",
	"Mean_OMP_Thres_2","Median_OMP_Thres_2",
	"Mean_OMP_Thres_4","Median_OMP_Thres_4")
colnames(results_df) <- Pheno_Names
rownames(results_df) <- algorithms



# ==================  MEAN ================ #
# ================== MEDIAN =============== #
mean_error <- median_error <- numeric(m)

# ==================  OMP ================= #
omp_error_2 <- omp_error_4 <- numeric(m)
# ================== SES ================== #
# == THreshold 0.01 
ses_error_001_2  <- ses_error_001_3  <- ses_error_001_4  <- ses_error_001_5  <- numeric(m)
# == Threshold 0.05 
ses_error_005_2  <- ses_error_005_3  <- ses_error_005_4  <- ses_error_005_5  <- numeric(m)

pheno_num =0


# =============================================================== #
# ========== FOR EVERY TRAIT OF INTEREST ======================== #
# =============================================================== #

for(pheno in list(Fat,Weight,Tag_Weight,Length_Width)){
pheno_num = pheno_num+1
print("pheno_num",pheno_num)
# =============================================================== #
# =========== FOR EVERY Observation that we leave out =========== #
# =============================================================== #
for (i in 1:dim(x)[1]){
    # Prepare dataset
    # == TRAIN
    Train_y    = pheno[-i]
    Train_X    = x[-i,]
    Train_id   = id[-i]
    
    # == TEST
    Test_y    = pheno[i]
    Test_X    = x[i,]
    Test_id   = id[i] 
  	
    # Impute missing Values with the most frequent value
    if (Mean_Impute=="y"){
        for (snp_mis in range(dim(Train_X)[2])){
            most_freq <- as.numeric(names(which.max(table(Train_X[snp_mis]))))
            Train_X[snp_mis][is.na(Train_X[snp_mis])] <- as.numeric(most_freq)
            Test_X[snp_mis][is.na(Test_X[snp_mis])]   <- as.numeric(most_freq)}
        Train_X = as.matrix(Train_X)
        Train_X <- apply(Train_X,2,as.numeric)
        Test_X = as.matrix(Test_X)
        Test_X <- apply(Test_X,2,as.numeric)}
    # ================================  MEAN =============================================== # 
	mod0 <- lmer(Train_y ~ 1 + (1|Train_id), REML = FALSE )   ## REML = TRUE
	b0 <- as.vector( coef(mod0)$Train_id )
	mean_error[i] <- sum( ( Test_y -b0 )^2 ) / m  
	# ================================= OMP ================================================ #
    # == Thres 2,4
    omp_thres = 2 
    Error = my_omp_fun2(Train_y,Train_X,Train_id,Test_y, Test_X, Test_id,omp_thres)
    omp_error_2[i] = Error[1]
    omp_error_4[i] = Error[2]
    #--------------------------------------------------------------------------------------- #

	# ================================= SES ================================================ #
    # ============== Threshold 0.01 =================== #
    ses_thres = 0.01 
    ses_k = 2
    Error = my_ses_fun(fam_num,Train_y,Train_X,Train_id,Test_y, Test_X, Test_id,ses_thres,ses_k)
    ses_error_001_2[i]  = Error
    ses_k = 3
    Error = my_ses_fun(fam_num,Train_y,Train_X,Train_id,Test_y, Test_X, Test_id,ses_thres,ses_k)
    ses_error_001_3[i]  = Error
    ses_k = 4
    Error = my_ses_fun(fam_num,Train_y,Train_X,Train_id,Test_y, Test_X, Test_id,ses_thres,ses_k)
    ses_error_001_4[i]  = Error
    ses_k = 5
    Error = my_ses_fun(fam_num,Train_y,Train_X,Train_id,Test_y, Test_X, Test_id,ses_thres,ses_k)
    ses_error_001_5[i]  = Error
    # ============== Threshold 0.05 =================== #
    ses_thres = 0.05 
    ses_k = 2
    Error = my_ses_fun(fam_num,Train_y,Train_X,Train_id,Test_y, Test_X, Test_id,ses_thres,ses_k)
    ses_error_005_2[i]  = Error
    ses_k = 3
    Error = my_ses_fun(fam_num,Train_y,Train_X,Train_id,Test_y, Test_X, Test_id,ses_thres,ses_k)
    ses_error_005_3[i]  = Error
    ses_k = 4
    Error = my_ses_fun(fam_num,Train_y,Train_X,Train_id,Test_y, Test_X, Test_id,ses_thres,ses_k)
    ses_error_005_4[i]  = Error
    ses_k = 5
    Error = my_ses_fun(fam_num,Train_y,Train_X,Train_id,Test_y, Test_X, Test_id,ses_thres,ses_k)
    ses_error_005_5[i]  = Error}
    # --------------------------------------------------------------------------------------- #  

results = c(mean(mean_error),median(mean_error),
	mean(ses_error_001_2),median(ses_error_001_2),
	mean(ses_error_001_3),median(ses_error_001_3),
	mean(ses_error_001_4),median(ses_error_001_4),
	mean(ses_error_001_5),median(ses_error_001_5),
	mean(ses_error_005_2),median(ses_error_005_2),
	mean(ses_error_005_3),median(ses_error_005_3),
	mean(ses_error_005_4),median(ses_error_005_4),
	mean(ses_error_005_5),median(ses_error_005_5),	
	mean(omp_error_2),median(omp_error_2),
	mean(omp_error_4),median(omp_error_4))
results_df[,pheno_num] = results
write.table(results, file=paste("Pheno_",pheno_num,"_RESULTS_CV.tsv",sep=""), quote=FALSE, sep='\t', col.names = NA)}

write.table(results_df, file="Cross_Validation_Tsamar/RESULTS_CV2.tsv", quote=FALSE, sep='\t', col.names = NA)

