######## cox lasso model for feature selection ###
#Load required libraries
library(caret)
library(ggfortify)
library(survival) 
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(MASS)
library(survivalROC)
library(glmnet)
library(ROCR) 
library(rms) 
library(Hmisc) 
args <- commandArgs(TRUE)

setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Primary_grade4_survival/")

set.seed(7)

data <- read.table("185_genes_data_with_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
data <- read.table("185_genes_minus11_genes_data_with_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#data <- read.table("185_genes_minus20_genes_data_with_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
#data <- read.table("185_genes_minus35_genes_data_with_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)


head(data[1:25],2)
dim(data)


head(data[25],2)

data1 <- subset(data, OS_time!="NA")
dim(data1)
data1 <- subset(data1, OS_month > 0)
dim(data1)
data1 <- subset(data1, Death_Status!="NA")
dim(data1)


clin <- data1[1:25]
head(clin,2)
Exp <- data1[26:ncol(data1)]
head(Exp[1:20],2)

Exp_t <- as.data.frame(t(Exp))



#write.table(cbind("set","Beta_Coeff","HR with 95% CI","C-Index","P-value","P-value1"), file="Train_Random11_genes_PI_results.txt",row.names=F,col.names=F,sep = '\t');
#write.table(cbind("set","Beta_Coeff","HR with 95% CI","C-Index","P-value","P-value1"), file="Ext_Test_Random11_genes_PI_results.txt",row.names=F,col.names=F,sep = '\t');

write.table(cbind("set","Beta_Coeff","HR with 95% CI","C-Index","P-value","P-value1"), file="Train_Random8_genes_PI_results.txt",row.names=F,col.names=F,sep = '\t');
write.table(cbind("set","Beta_Coeff","HR with 95% CI","C-Index","P-value","P-value1"), file="Ext_Test_Random8_genes_PI_results.txt",row.names=F,col.names=F,sep = '\t');



#for(i in 1:100) {

  for(n in 1:100) {
    
print(n)
  
#Exp1_t <- Exp_t[sample(nrow(Exp_t), 11), ]
Exp1_t <- Exp_t[sample(nrow(Exp_t), 8), ]
Exp1_train <- as.data.frame(t(Exp1_t))
#Exp1 <- Exp[sample(ncol(Exp), 3), ]
head(Exp1 )


sel_ftrs <- as.data.frame(names(Exp1_train))
names(sel_ftrs) <- c("ID" )
sel_ftrs

######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
train_feature_mat<- cbind(data1["OS_month"],data1["Death_Status"],Exp1_train)
head(train_feature_mat,2)
dim(train_feature_mat)

############ remove where OS.time=NA ############
train_feature_mat1<-subset(train_feature_mat,OS_month!="NA")
train_feature_mat1<-subset(train_feature_mat1,OS_month!=0)


#train_feature_mat2<- cbind(data2["OS_month"],data2["Group"],sel_train)
# save files with selected genes & survival information #########
#write.table(cbind("ID"=rownames(train_feature_mat1), train_feature_mat1),file="sel_train.txt",sep="\t",quote=F, row.names=F)


#Create prognostic index for pathway
tr <- train_feature_mat1[3:ncol(train_feature_mat1)]

E= length(tr)
E

PI_tr = 0
for(i in seq(from=1, to=E ,by=1))
{
  #PI_tr= PI_tr+((tr[,i])*(sel_features_results[i,1]))
  PI_tr= PI_tr+((tr[,i])* 1)
}



train_feature_mat1$PI<-PI_tr
head(train_feature_mat1)

#save selected data with PI value
#write.table(cbind("ID"=rownames(LUSC_C_tr), LUSC_C_tr),file="train_11_genes_with_PI.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(LUSC_C_tr), LUSC_C_tr),file="train_8_genes_with_PI.txt",sep="\t",quote=F, row.names=F)

tr_PI <- as.data.frame(train_feature_mat1$PI)
rownames(tr_PI) <- rownames(train_feature_mat1)
colnames(tr_PI) <- c("PI")
#write.table(cbind("ID"=rownames(tr_PI ), tr_PI ),file="train_11genes_based_PI.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(tr_PI ), tr_PI ),file="train_8genes_based_PI.txt",sep="\t",quote=F, row.names=F)


######################################## Survival Object ############################################
surv_object_tr <- Surv(time = train_feature_mat1$OS_month, event = train_feature_mat1$Death_Status)
surv_object_tr
dim(surv_object_tr )


dim(tr)
head(tr_PI)
mean(tr_PI$PI)

# ###survival analysis: fits cox ph model to find HR for PI
fit_tr <- survfit(surv_object_tr~(tr_PI$PI>mean(tr_PI$PI)), data=tr_PI)
summary(fit_tr)

#fit.coxph_tr <- coxph(surv_object_tr ~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
fit.coxph_tr <- coxph(surv_object_tr ~(tr_PI$PI>mean(tr_PI$PI)), data=tr_PI)
summary(fit.coxph_tr)

tr_res <- summary(fit.coxph_tr)
coeff <- round(tr_res$coefficients[1],2)
HR <- round(tr_res$conf.int[1],2)
int1 <- round(tr_res$conf.int[3],2)
int2 <- round(tr_res$conf.int[4],2)
CI <- round (tr_res$concordance[1],2)
pval <- tr_res$sctest[3]
pval1 <- format(pval, scientific=T)
pval1

HR1 <- paste0(HR , " (",  int1,  " - ", int2, ") " )
HR1

tr_res1 <- cbind(coeff, HR1, CI, pval,pval1)
row.names(tr_res1) <- paste0("Random_set_",n)
#row.names(tr_res1) <- c("set")
tr_res1

#write.table(tr_res1,file="Train_Random11_genes_PI_results.txt",row.names=T,col.names=F,sep = '\t',append = T);#output file
write.table(tr_res1,file="Train_Random8_genes_PI_results.txt",row.names=T,col.names=F,sep = '\t',append = T);#output file


### External  Test Data ##########



#################### External Validation dataset ############
Test_data <- read.table("Ext_Primary_GBM_both_trt_Clin_data_lts_mts_sts_with_185_QN_data.txt", sep="\t", header=T, check.names = F, row.names = 1)
dim(Test_data)
head(Test_data[1:30],2)



#data preparation with selected features 
#sel_train <- as.data.frame(data1[,colnames(data1) %in% c(row.names(features)), ])
#sel_train <- as.data.frame(data1[,colnames(data1) %in% c(sel_features_results$ID), ])
sel_test <- as.data.frame(Test_data[,colnames(Test_data) %in% c(sel_ftrs$ID), ])
head(sel_test ,2)
dim(sel_test )

######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
test_feature_mat<- cbind(Test_data["OS_month"],Test_data["Death_event"],sel_test )
head(test_feature_mat,2)
dim(test_feature_mat)

############ remove where OS.time=NA ############
test_feature_mat1<-subset(test_feature_mat,OS_month!="NA")
test_feature_mat1<-subset(test_feature_mat1,OS_month!=0)
dim(test_feature_mat1)

#Create prognostic index for pathway
te3 <- test_feature_mat1[3:ncol(test_feature_mat1)]

E3= length(te3)
E3


PI_te3=0
for(i in seq(from=1, to=E3 ,by=1))
{
  #PI_tr3= PI_tr3+((tr3[,i])*(sel_features_results[i,1]))
  PI_te3= PI_te3+((te3[,i])*1)
}

test_feature_mat1$PI<-PI_te3
head(test_feature_mat1)

#save selected data with PI value
#write.table(cbind("ID"=rownames(LUSC_C_tr3), LUSC_C_tr3),file="Ext_test_LTS_STS_with_11_genes_PI.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(LUSC_C_tr3), LUSC_C_tr3),file="Ext_test_LTS_MTS_STS_with_8genes_PI.txt",sep="\t",quote=F, row.names=F)


te_PI3 <- as.data.frame(test_feature_mat1$PI)
rownames(te_PI3) <- rownames(test_feature_mat1)
colnames(te_PI3) <- c("PI")
head(te_PI3)
#write.table(cbind("ID"=rownames(tr_PI3 ), tr_PI3 ),file="te_PI3",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(tr_PI3 ), tr_PI3 ),file="test_PI_based_on_8genes.txt",sep="\t",quote=F, row.names=F)

######################################## Survival Object ############################################
surv_object_te<- Surv(time = test_feature_mat1$OS_month, event = test_feature_mat1$Death_event)
surv_object_te
dim(surv_object_te )


# ###survival analysis: fits cox ph model to find HR for PI
fit_te <- survfit(surv_object_te~(te_PI3$PI>mean(te_PI3$PI)), data=te_PI3)
#fit_tr <- survfit(surv_object_tr~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
summary(fit_te)

#fit.coxph_tr <- coxph(surv_object_tr ~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
fit.coxph_te <- coxph(surv_object_te ~(te_PI3$PI>mean(te_PI3$PI)), data=te_PI3)
summary(fit.coxph_te)


te_res <- summary(fit.coxph_te)
coeff_te <- round(te_res$coefficients[1],2)
HR_te <- round(te_res$conf.int[1],2)
int1_3 <- round(te_res$conf.int[3],2)
int2_3 <- round(te_res$conf.int[4],2)
CI_te <- round (te_res$concordance[1],2)
pval_te <- te_res$sctest[3]
pval1_te <- format(pval_te, scientific=T)


HR1_te <- paste0(HR_te , " (",  int1_3,  " - ", int2_3, ") " )


te_res1 <- cbind(coeff_te, HR1_te, CI_te, pval_te , pval1_te)
te_res1

#row.names(te_res1) <- c("set")
row.names(te_res1) <- paste0("Random_set_",n)
te_res1


#write.table(te_res1,file="Ext_Test_Random11_genes_PI_results.txt",row.names=T,col.names=F,sep = '\t',append = T);#output file
write.table(te_res1,file="Ext_Test_Random8_genes_PI_results.txt",row.names=T,col.names=F,sep = '\t',append = T);#output file



}

