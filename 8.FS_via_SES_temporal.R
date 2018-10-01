library(MXM)
library(lme4)
library(grid)
library(gridExtra)
library(gtable)
library(ggplot2)
library(webshot)
library(cluster)



# ============= MODES ============== #
With_impute = "y"
thres = 0.01
thres_n = "001"
condition_set = 2
# ================================== #


# ============================ LOAD DATA =================================== #
mydata = read.csv("Dataset.tsv",header = TRUE,sep="\t")
Image_Dir = paste("Feature_Selection/Ses_thres_",thres_n,"_k_",condition_set,"/",name_file,"/",sep="")
# ========================================================================= #



attach(mydata)
Pheno    = as.data.frame(mydata[c(3,6,7,4)])#Phenotypes

All_Features = as.data.frame(mydata[8:dim(mydata)[2]])
Names = colnames(Pheno)
IDS = paste("Family",Fam_ID)
group = Fam_ID

counter = 0
dataset = as.matrix(All_Features)
Features  = All_Features
Family_Id = IDS



for (pheno in Pheno) {
    counter = counter +1
    name = Names[counter]
    target = pheno
    
    m1 = SES.temporal(target, reps = NULL, group, dataset, max_k = condition_set, threshold = thres, 
                          test = "testIndLMM",ncores = 1)
    
    Sign     = m1@signatures
    Summ     = capture.output(summary(m1))
    N_Sign   = nrow(Sign)
    Var_l    = length(colnames(Sign))
    All_vars1 = colnames(Features)[unique(c(Sign))]
  

    final_vars2 = colnames(All_Features)[m1@selectedVars]
    target = c(target)
    form = paste("target ~ ", paste(final_vars2,collapse="+"), "+ (1|group)")
    dataset_lmer = as.data.frame(dataset)
    m2 = lmer(form,data=dataset_lmer,REML=FALSE)

    d<-data.frame(pheno,pred = fitted(m2),resid = resid(m2))
    f_bic = round(BIC(m2), digits = 2)

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
    

    All_vars = All_vars1
    # Create Tables with the different signatures produced by SES
    for (i in 1:N_Sign){
        var_indx = as.vector(Sign[i,])
        vars     = colnames(Features)[var_indx]
        var_t = vector()
        for (var_i in 1:length(vars)){
        	var_t[var_i] = paste(unlist(strsplit(vars[var_i], "[.]"))[1:2],collapse=":")
        } 
        pvals    = exp(m1@pvalues[var_indx])
        # ======== PLOT MATRIX =================== #
        Fat_T = matrix(data = c(var_t,pvals),ncol=2)
        colnames(Fat_T) = c("Variables" , "P-value")
        write.csv(Fat_T, paste("Feature_Selection/",With_impute,thres_n,name,i,".tsv",sep=""), sep='\t')
        t1 <- tableGrob(Fat_T)
        title <- textGrob(name,gp=gpar(fontsize=30))
        padding <- unit(5,"mm")
        table <- gtable_add_rows( t1, heights = grobHeight(title) + padding, pos = 0)
        table <- gtable_add_grob( table, title, 1, 1, 1, ncol(table))
        grid.newpage()
        Basi =paste(Image_Dir,name,"/",i,"_SESTEMP_",name, sep = "")
        print(Basi)
        png(paste(Basi,".png", sep = ""))
        grid.draw(table)
        dev.off()
    }
    # Create a boxplot and a violin plot for all selected variables
    for (variable in All_vars) {
        Ref = unlist(strsplit(variable,"[.]"))[4]
        Alt = unlist(strsplit(variable,"[.]"))[5]
        HomR= paste(Ref,"/",Ref,sep="")
        Hetr= paste(Ref,"/",Alt,sep="")
        HomA= paste(Alt,"/",Alt,sep="")


        Title = paste(name," vs ",variable)

        plot_title =  paste(unlist(strsplit(variable, "[.]"))[1:2],collapse=":")

        v1 = factor(mydata[,variable], levels = c('0','1','2'),ordered = TRUE)
        ggplot(mydata, aes(v1,pheno)) +
                    geom_violin() +
                    geom_point(aes(colour=Family_Id))+
                    ggtitle(plot_title) +
                    scale_x_discrete(labels=c("0" = HomR, "1" = Hetr,
                              "2" = HomA)) +
                    theme(text = element_text(size=20),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5)) +
                    ylab(name)
            
        Basic_T = paste(Image_Dir,name,'/',Title, sep = "")
        ggsave(paste(Basic_T,'_violin.pdf', sep = ""))

        ggplot(mydata, aes(v1,pheno)) +
                geom_boxplot() +
                geom_jitter(aes(colour=Family_Id))+
                ggtitle(plot_title) +
                scale_x_discrete(labels=c("0" = HomR, "1" = Hetr,
                              "2" = HomA)) +
                theme(text = element_text(size=20),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5)) +
                ylab(name)
        ggsave(paste(Basic_T,'.pdf', sep = ""))


  }
}
