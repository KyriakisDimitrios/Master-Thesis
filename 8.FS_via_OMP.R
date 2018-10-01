library(lme4)
library(grid)
library(gridExtra)
library(gtable)
library(ggplot2)
library(webshot)
library(plotly)
library(ggfortify)
library(cluster)





# ============= MODES ============== #
thres = 4
# ================================== #


# ============================ LOAD DATA =================================== #
mydata = read.csv("Dataset.tsv",header = TRUE,sep="\t")
mydata[ mydata == "-" ] <- NA

Image_Dir = paste("Feature_Selection/OMP_thres_",thres,"/",sep="")
# ========================================================================= #




Pheno    = as.data.frame(mydata[c(3,6,7)])#Phenotypes
attach(Pheno)
All_Features = as.data.frame(mydata[8:dim(mydata)[2]])# SNPs
Names = colnames(Pheno)
IDS_S = paste("Family",mydata$Fam_ID)
IDS = mydata$Fam_ID
counter = 0
group = mydata$Fam_ID
dataset = as.matrix(All_Features)



#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------


# =================================== MODELS =====================================
myfun<- function(i){
  return(cor(resid(m1),i))
}
#---------------------------------------------------------------------------------

col_names = colnames(All_Features)
pheno_num =1
for (pheno_num in 1:4){
    name  = Names[pheno_num] 
    print(name)
    pheno = Pheno[,pheno_num]
    
    # Initialize
    form = "pheno ~+ (1|IDS)"  
    m1 = lmer(form,REML=FALSE)
    BIC_old = BIC(m1)
    Diffe = 5
    counter= 0
    list_var= c()
    
    #  OMP Iteration
    while (Diffe >thres & BIC_old >0 )
    {
      counter = counter+1
      cor_vals = sapply(All_Features,myfun)
      a = which.max(cor_vals)
      list_var[counter]  = col_names[a]
      form = paste("pheno ~ ", paste(list_var,collapse="+"), "+ (1|IDS)")
      m1 = lmer(form,data=All_Features,REML=FALSE)
      BIC_new = BIC(m1)
      Diffe = BIC_old - BIC_new
      BIC_old = BIC_new
      
    }
    
    # Final model
    final_vars = list_var[-length(list_var)]

    
    form = paste("pheno ~ ", paste(final_vars,collapse="+"), "+ (1|IDS)")
    m4 = lmer(form,data=All_Features,REML=FALSE)


    d<-data.frame(pheno,pred = fitted(m4),resid = resid(m4))

    # ================= FITTED VS OBESSRVED
    f_bic = round(BIC(m4), digits = 2)
    
    # ==================== FITTED VS RESIDUALS
    ggplot(d, aes(x = pred, y = resid)) +
          geom_segment(aes(xend =pred, yend = resid), alpha = .2) +
          geom_point(aes(y = resid), shape = 1,color="blue") +
          geom_hline(yintercept=0, linetype="dashed", color = "black")+
          ggtitle(paste("Fitted vs Residuals\n","BIC = ",f_bic,sep=""))+
          xlab("Fitted") + ylab("Residuals (type='pearson')")+
          theme(plot.title = element_text(face="bold",hjust=0.5),
                axis.title.x = element_text( size=14),
            axis.title.y = element_text( size=14))
    plot_name =paste(Image_Dir,name,'_Fitted_vs_Errors.jpg',sep="")
    ggsave(plot_name)


    # ==================== FITTED VS OBSERVED
    ggplot(d, aes(x = pheno, y = pred)) +
          geom_smooth(method = "lm", se = FALSE, color = "lightgrey") +
          geom_segment(aes(xend =pheno, yend = pred), alpha = .2) +
          geom_point(aes(y = pred), shape = 1) +
          ggtitle("Observed vs Predicted")+
          xlab("Observed") + ylab("Predicted")+
          geom_abline() +
          theme(plot.title = element_text(face="bold",hjust=0.5),
                axis.title.x = element_text( size=14),
            axis.title.y = element_text( size=14))
    plot_name =paste(Image_Dir,name,'_Fitted_vs_Real.png',sep="") 
    ggsave(plot_name)
    
    #=================================================================
    # Continue if no vars
    if (length(final_vars)==0) next

    vars     = final_vars
    
    # # ======== PLOT MATRIX With SElected SNPs =================== #
    Fat_T = matrix(data = c(vars),ncol=1)
    colnames(Fat_T) = c("Variables")
    t1 <- tableGrob(Fat_T)
    title <- textGrob(name,gp=gpar(fontsize=30))
    padding <- unit(5,"mm")
    table <- gtable_add_rows( t1, heights = grobHeight(title) + padding, pos = 0)
    table <- gtable_add_grob( table, title, 1, 1, 1, ncol(table))
    grid.newpage()
    Basi =paste(Image_Dir,name,"/OMP_",name, sep = "")
    png(paste(Basi,".png", sep = ""))
    grid.draw(table)
    dev.off()
    
    Family_Id = IDS_S
    # Boxplot and Violin plot for all selected SNPs
    for (variable in final_vars) {
        Title = paste(name," vs ",variable)
        v1 = factor(mydata[,variable])
        
        ggplotly(ggplot(mydata, aes(v1,pheno)) +
                geom_violin() +
                geom_point(aes(colour=Family_Id))+
                xlab(variable) +
                theme(text = element_text(size=20)) +
                ylab(name))
        
        Basic_T = paste(Image_Dir,name,'/',Title, sep = "")
        ggsave(paste(Basic_T,'_violin.png', sep = ""))

        ggplotly(ggplot(mydata, aes(v1,pheno)) +
            geom_boxplot() +
            geom_jitter(aes(colour=IDS_S))+
            xlab(variable) +
            theme(text = element_text(size=20)) +
            ylab(name))
        ggsave(paste(Basic_T,'.png', sep = ""))
    }

}

