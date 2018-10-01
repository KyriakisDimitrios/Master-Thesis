# ========= LIBRARIES =========== #
library(MXM)
library(lme4)


# ======== FUNCTIONS =========== # 

# ====================== OMP ============================== #
# ====================== OMP ============================== #
my_omp_fun <- function(Train_y,Train_X,Train_id,Test_y, Test_X, Test_id,omp_thres){
    # Initialize
    m1 = lmer(Train_y~1+(1|Train_id),REML=FALSE)
    BIC_old = BIC(m1)
    Diffe = 5
    counter= 0
    list_var= c()
    Train_X = as.data.frame(Train_X)
    #  OMP Iteration
    while (Diffe >omp_thres)
    {
      col_names = colnames(Train_X)
      counter = counter+1
      cor_vals = c(cor(resid(m1),Train_X))
      a = which.max(cor_vals)
      list_var[counter]  = col_names[a]
      form = paste("Train_y ~ ", paste(list_var,collapse="+"), "+ (1|Train_id)")
      m1 = lmer(form,data=Train_X,REML=FALSE)
      BIC_new = BIC(m1)
      Diffe = BIC_old - BIC_new
      BIC_old = BIC_new
      
    }
    
    # Final model
    sel = list_var[-length(list_var)]
    form = paste("Train_y ~ ", paste(sel,collapse="+"), "+ (1|Train_id)")
    mod1 = lmer(form,data=Train_X,REML=FALSE)
    b1 <- as.matrix( coef(mod1)$Train_id )
    Error <- sum( ( Test_y - c(1, Test_X[sel]) %*% t(b1) )^2 ) / m
    return(Error)
}

# ====================== SES ============================== #
my_ses_fun <- function(fam_num,Train_y,Train_X,Train_id,Test_y, Test_X, Test_id,ses_thres,ses_k){
    #  SES
	sel <- SES.temporal(target = Train_y, dataset = Train_X, group = Train_id, test= "testIndLMM",max_k = ses_k,threshold = ses_thres)@selectedVars
	print(length(sel))
    mod1 <- lmer(Train_y ~ Train_X[ , sel ] + (1|Train_id), REML = FALSE )  ## REML = TRUE 
	b1 <- as.matrix( coef(mod1)$Train_id )
	Error <- sum( ( Test_y - c(1, Test_X[sel]) %*% t(b1) )^2 ) / m
    print(Error)    
    return(Error)
}
# ------------------------------------------------------------------------------------------- #

# ====================== OMP ============================== #
my_omp_fun2 <- function(Train_y,Train_X,Train_id,Test_y, Test_X, Test_id,omp_thres){
    # Initialize
    m1 = lmer(Train_y~1+(1|Train_id),REML=FALSE)
    BIC_old = BIC(m1)
    Diffe = 5
    counter= 0
    list_var= c()
    Train_X = as.data.frame(Train_X)
    #  OMP Iteration
    Diffe_list=c()
    while (Diffe >=omp_thres)
    {
        col_names = colnames(Train_X)
        counter = counter+1
        cor_vals = c(cor(resid(m1),Train_X))
        a = which.max(cor_vals)
        list_var[counter]  = col_names[a]
        form = paste("Train_y ~ ", paste(list_var,collapse="+"), "+ (1|Train_id)")
        m1 = lmer(form,data=Train_X,REML=FALSE)
        BIC_new = BIC(m1)
        Diffe = BIC_old - BIC_new
        BIC_old = BIC_new
        Diffe_list =c(Diffe_list,Diffe) 
    }
    
    # Final model
    sel_bic4 = list_var[1:min(which(Diffe_list<4))]
    sel_bic4 = sel_bic4[-length(sel_bic4)]
    sel_bic2 = list_var[-length(list_var)]
    
    # OMP 2
    if (length(sel_bic2)==0){
        mod0 <- lmer(Train_y ~ 1 + (1|Train_id), REML = FALSE )   ## REML = TRUE
        b0 <- as.vector( coef(mod0)$Train_id )
        Error_omp2 <- sum( ( Test_y -b0 )^2 ) / m  
    }else{
        form = paste("Train_y ~ ", paste(sel_bic2,collapse="+"), "+ (1|Train_id)")
        mod1 = lmer(form,data=Train_X,REML=FALSE)
        # mod1 <- lmer(Train_y ~ Train_X[ , sel ] + (1|Train_id), REML = FALSE )  ## REML = TRUE 
        b1 <- as.matrix( coef(mod1)$Train_id )
        Error_omp2 <- sum( ( Test_y - c(1, Test_X[sel_bic2]) %*% t(b1) )^2 ) / m}
        
    # OMP 4
    if (length(sel_bic4)==0){
        mod0 <- lmer(Train_y ~ 1 + (1|Train_id), REML = FALSE )   ## REML = TRUE
        b0 <- as.vector( coef(mod0)$Train_id )
        Error_omp4 <- sum( ( Test_y -b0 )^2 ) / m  
    }else{
        form = paste("Train_y ~ ", paste(sel_bic4,collapse="+"), "+ (1|Train_id)")
        mod1 = lmer(form,data=Train_X,REML=FALSE)
        # mod1 <- lmer(Train_y ~ Train_X[ , sel ] + (1|Train_id), REML = FALSE )  ## REML = TRUE 
        b1 <- as.matrix( coef(mod1)$Train_id )
        Error_omp4 <- sum( ( Test_y - c(1, Test_X[sel_bic4]) %*% t(b1) )^2 ) / m }
    Error = c(Error_omp2,Error_omp4)
    return(Error)
}