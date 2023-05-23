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
#data <- read.table("185_genes_minus11_genes_data_with_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

head(data[1:25],2)
dim(data)


head(data[25],2)

data1 <- subset(data, OS_time!="NA")
dim(data1)
data1 <- subset(data1, OS_month > 0)
dim(data1)
data1 <- subset(data1, Death_Status!="NA")
dim(data1)

#data1 <- subset(data1, OS < 0)

# Select features using Lasso Regression 
surv_object2 <- Surv(time = data1$OS_month, event = data1$Death_Status)
set.seed(10)
cvfit1 <- cv.glmnet(as.matrix(data1[26:ncol(data1)]), 
                    surv_object2, # create survival object from the data
                    family = "cox", # specify Cox PH model
                    type.measure = "C", 
                    nfolds = 5, 
                    alpha = 1, # lasso: alpha = 1; ridge: alpha=0
                    maxit = 1000)

lambda_min <- cvfit1$lambda.min
lambda_min

plot(cvfit1)
cvfit1

#plot(cvfit1)

jpeg(file="Cox_Lasso_Regression_lamda_plot.jpeg", units="in", width=10, height=10, res=350)
plot(cvfit1)
dev.off()

# plot lambda coefficinets
plot(cvfit1, xvar = "lambda", yvar= "Measure", label = TRUE)

plot(cvfit1, xvar = "dev", label = TRUE)

est.coef = coef(cvfit1, s = cvfit1$lambda.min) # returns the p length coefficient vector
head(est.coef)
# of the solution corresponding to lambda 
active.k = which(est.coef != 0)
active.k

#Extract coefficinet values
active.k.vals = est.coef[active.k]
active.k.vals

key_variables <- as.data.frame(est.coef[est.coef[,1]!=0,])
colnames(key_variables) <- c("coeff")
key_variables <- round(key_variables,3)
key_variables 
dim(key_variables)

#write into a file
#write.table(cbind("ID"=rownames(key_variables), key_variables),file="Lasso_8key_variables_based_on_STS_MTS_LTS_data_all_responsive_genes.txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(key_variables), key_variables),file="Lasso_11key_variables_based_on_STS_MTS_LTS_data_all_responsive_genes.txt",sep="\t",quote=F, row.names=F)


############## Develop Prognostic Index  Model based on selected features #######

######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
#sel_features_results<-read.table("Lasso_8key_variables_based_on_STS_MTS_LTS_data_all_responsive_genes.txt",header =TRUE, sep = "\t", row.names=1, check.names = FALSE)

sel_features_results<-read.table("Lasso_11key_variables_based_on_STS_MTS_LTS_data_all_responsive_genes.txt",header =TRUE, sep = "\t", row.names=1, check.names = FALSE)

head(sel_features_results)
dim(sel_features_results)


#training data preparation with selected features 
#sel_train <- as.data.frame(data1[,colnames(data1) %in% c(sel_features_results$ID), ])
sel_train <- as.data.frame(data1[,colnames(data1) %in% c(row.names(sel_features_results)), ])
head(sel_train,2)
dim(sel_train)

######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
train_feature_mat<- cbind(data1["OS_month"],data1["Death_Status"],sel_train)
head(train_feature_mat,2)
dim(train_feature_mat)

############ remove where OS.time=NA ############
train_feature_mat1<-subset(train_feature_mat,OS_month!="NA")
train_feature_mat1<-subset(train_feature_mat1,OS_month!=0)

# save files with selected genes & survival information #########
write.table(cbind("ID"=rownames(train_feature_mat1), train_feature_mat1),file="sel_train_11_features.txt",sep="\t",quote=F, row.names=F)


#Create prognostic index 
LUSC_C_tr=train_feature_mat1
tr <- LUSC_C_tr[3:ncol(LUSC_C_tr)]

E= length(tr)
E

head(tr,4)
tr[,2]
sel_features_results[2,1]
tail(sel_features_results,2)

PI_tr = 0
for(i in seq(from=1, to=E ,by=1))
{
  PI_tr= PI_tr+((tr[,i])*(sel_features_results[i,1]))
  #PI_tr= PI_tr+((tr[,i])* 1)
}



# add PI as new column to the data
LUSC_C_tr$PI<-PI_tr
head(LUSC_C_tr)

#save selected data with PI value
write.table(cbind("ID"=rownames(LUSC_C_tr), LUSC_C_tr),file="train_11_genes_with_PI.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(LUSC_C_tr), LUSC_C_tr),file="train_8_genes_with_PI.txt",sep="\t",quote=F, row.names=F)

tr_PI <- as.data.frame(LUSC_C_tr$PI)
rownames(tr_PI) <- rownames(LUSC_C_tr)
colnames(tr_PI) <- c("PI")
write.table(cbind("ID"=rownames(tr_PI ), tr_PI ),file="train_11genes_based_PI.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(tr_PI ), tr_PI ),file="train_8genes_based_PI.txt",sep="\t",quote=F, row.names=F)


######################################## Survival Object ############################################
surv_object_tr <- Surv(time = LUSC_C_tr$OS_month, event = LUSC_C_tr$Death_Status)
surv_object_tr
dim(surv_object_tr )


# Fit survival model for KM plots
fit_tr <- survfit(surv_object_tr~(tr_PI$PI>mean(tr_PI$PI)), data=tr_PI)
summary(fit_tr)

##### survival analysis: fits cox ph model to find HR for PI
fit.coxph_tr <- coxph(surv_object_tr ~(tr_PI$PI>mean(tr_PI$PI)), data=tr_PI)
summary(fit.coxph_tr)

#extract important survival results information from the results object
tr_res <- summary(fit.coxph_tr) #get summary
coeff <- round(tr_res$coefficients[1],2) #beta coeff
HR <- round(tr_res$conf.int[1],2) #hazard ratio
int1 <- round(tr_res$conf.int[3],2) #HR first interval 
int2 <- round(tr_res$conf.int[4],2) #HR 2nd interval 
CI <- round (tr_res$concordance[1],2) # concornace index
pval <- tr_res$sctest[3] #pvalue
pval1 <- format(pval, scientific=T) #p-value in scientifc notaion


HR1 <- paste0(HR , " (",  int1,  " - ", int2, ") " ) #create HR with 95% CI

tr_res1 <- cbind(coeff, HR1, CI, pval1) #create table of results
names(tr_res1) <- c("set","Beta_Coeff","HR with 95% CI","C-Index","P-value")
row.names(tr_res1) <- c("Training_data")

#save KM survival plot
#jpeg(file="KM_plot.jpeg", units="in", width=10, height=10, res=300)
pp <-  ggsurvplot(fit_tr, data=tr_PI,
                  #pval=TRUE,
                  risk.table=TRUE,
                  tables.height = 0.3, #add risk table & height
                  xlab="Time in Months",
                  risk.table.col="strata", break.time.by = 12,
                  conf.int = F, censor = TRUE,
                  #title= paste0("KM plot based on PI/Risk score"),
                  surv.median.line = "hv", # Add medians survival
                  palette=c("dodgerblue2", "red"), #add desired color
                  size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                  #legend.title = paste0(pathway),
                  legend.labs = c("Less than mean " , " more than Mean"),
                  font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
pp

# customised the plot: add HR, CI and P-value on the KM plot
pp$plot <- pp$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 1, hjust = 1,
    #label = "HR = 0.9 \n p < 0.001",
    #label =  paste0("HR = ", HR1, \n, "p-val <", pval),
    label =  paste0("HR = ", HR1, "\n",  "p-val = ", pval1,  "\n",  "C-Index = ", CI),
    size = 5)

# now plot
pp


jpeg("KM_plot_11_genes_with_coeff.jpg", units="in", width=10, height=10, res=350)
#jpeg("KM_plot_8_genes_with_coeff.jpg", units="in", width=10, height=10, res=350)
print(pp, newpage = FALSE)
dev.off()

########## Draw KM plot Based on Quantile ranges of PI ######

fit_tr_Q <- survfit(surv_object_tr~ cut(tr_PI$PI, quantile(tr_PI$PI)), data=tr_PI)
summary(fit_tr_Q)


fit.coxph_tr_Q <- coxph(surv_object_tr ~ cut(tr_PI$PI, quantile(tr_PI$PI)), data=tr_PI)
summary(fit.coxph_tr_Q)

pp_Q <- ggsurvplot(fit_tr_Q, data=tr_PI,
                   #pval=TRUE,
                   risk.table=TRUE,
                   xlab="Time in Months",
                   tables.height = 0.35,
                   surv.median.line = "hv", # Add medians survival
                   palette=c("green","dodgerblue2", "pink","red"), #add desired color
                   size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                   #title= paste0("KM plot based on Quantiles of PI/Risk score"),
                   legend.title = paste0("Quantiles of PI score"),
                   #legend.title = paste0(pathway),
                   #legend.labs = c("Q1" , "Q2", "Q3", "Q4"),
                   risk.table.col="strata", break.time.by = 6)
#fit_tr <- survfit(surv_object_tr~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)

pp_Q


tr_res_Q <- summary(fit.coxph_tr_Q)
coeff_Q <- round(tr_res_Q$coefficients[1],2)
HR_Q <- round(tr_res_Q$conf.int[1],2)
int1_Q <- round(tr_res_Q$conf.int[3],2)
int2_Q <- round(tr_res_Q$conf.int[4],2)
CI_Q <- round (tr_res_Q$concordance[1],2)
pval_Q <- tr_res_Q$sctest[3]
pval1_Q <- format(pval_Q, scientific=T)
pval1_Q

HR1_Q <- paste0(HR_Q , " (",  int2_Q,  " - ", int1_Q, ") " )
HR1_Q

tr_res1_Q <- cbind(coeff_Q, HR1_Q, CI_Q, pval1_Q)
tr_res1_Q
#save results as a file
write.table(cbind("ID"=rownames(tr_res1_Q), tr_res1_Q),file="tr_results_11genes_based_PI_based_on Q.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(tr_res1), tr_res1),file="3genes_tr_res_PI.txt",sep="\t",quote=F, row.names=F)

pp_Q$plot <- pp_Q$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 1, hjust = 1,
    #label = "HR = 0.9 \n p < 0.001",
    #label =  paste0("HR = ", HR1, \n, "p-val <", pval),
    label =  paste0("HR = ", HR1_Q, "\n",  "p-val = ", pval1_Q,  "\n",  "C-Index = ", CI_Q),
    size = 5)

# now plot
pp_Q

jpeg("KM_plot_train_11_genes_with_coeff_with_Quantiles.jpg", units="in", width=10, height=10, res=350)
#jpeg("KM_plot_train_8_genes_with_coeff_with_Quantiles.jpg", units="in", width=10, height=10, res=350)
print(pp_Q, newpage = FALSE)
dev.off()

##Create ROC plot for 1-,3- and 5-years survival time prediction
tr_roc1 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                       status       = LUSC_C_tr$Death_Status,
                       marker       = LUSC_C_tr$PI,
                       predict.time = 12,
                       method       = "KM", 
                       #lambda = lambda_min , 
                       span = NULL,
                       window ="symmetric")
tr_roc1

tr_roc3 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                       status       = LUSC_C_tr$Death_Status,
                       marker       = LUSC_C_tr$PI,
                       predict.time = 36,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc3

tr_roc4 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                       status       = LUSC_C_tr$Death_Status,
                       marker       = LUSC_C_tr$PI,
                       predict.time = 48,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc4
tr_roc5 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                       status       = LUSC_C_tr$Death_Status,
                       marker       = LUSC_C_tr$PI,
                       predict.time = 60,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc5
tr_roc6 <- survivalROC(Stime        = LUSC_C_tr$OS_month,
                       status       = LUSC_C_tr$Death_Status,
                       marker       = LUSC_C_tr$PI,
                       predict.time = 72,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")

tr_roc6


plot(tr_roc1$FP, tr_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red",
     ylab="TP",main= "AUC Curve for Survival Prediction")
lines(tr_roc3$FP, tr_roc3$TP, type="l", lty=2, col="blue")
lines(tr_roc5$FP, tr_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("3 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)

#save plot ######
jpeg(file="ROC_train_11_genes_with_coeff.jpeg", units="in", width=10, height=10, res=300)
#jpeg(file="ROC_train_8_genes.jpeg", units="in", width=10, height=10, res=350)
plot(tr_roc1$FP, tr_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red",
     ylab="TP",main= "AUC Curve for Survival Prediction")
lines(tr_roc3$FP, tr_roc3$TP, type="l", lty=2, col="blue")
lines(tr_roc5$FP, tr_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("3 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)
dev.off()

#################### External Validation dataset ############
Test_data <- read.table("Ext_Primary_GBM_both_trt_Clin_data_lts_mts_sts_with_185_QN_data.txt", sep="\t", header=T, check.names = F, row.names = 1)
dim(Test_data)
head(Test_data[1:30],2)



#Test data preparation with selected features 
sel_test <- as.data.frame(Test_data[,colnames(Test_data) %in% c(row.names(sel_features_results)), ])
head(sel_test,2)
dim(sel_test)

######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
test_feature_mat2<- cbind(data3["OS_month"],data3["Death_event"],sel_test3)
head(test_feature_mat2,2)
dim(test_feature_mat2)

############ remove where OS.time=NA ############
test_feature_mat3<-subset(test_feature_mat2,OS_month!="NA")
test_feature_mat3<-subset(test_feature_mat3,OS_month!=0)
dim(test_feature_mat3)

#train_feature_mat2<- cbind(data2["OS_month"],data2["Group"],sel_train)
# save files with selected genes & survival information #########
#write.table(cbind("ID"=rownames(train_feature_mat1), test_feature_mat3),file="sel_test.txt",sep="\t",quote=F, row.names=F)

#Create prognostic index for test data
LUSC_C_te3=test_feature_mat3
te3 <- LUSC_C_te3[3:ncol(LUSC_C_te3)]

E3= length(te3)
E3


PI_te3=0
for(i in seq(from=1, to=E3 ,by=1))
{
  PI_te3= PI_te3+((te3[,i])*(sel_features_results[i,1]))
  #PI_tr3= PI_tr3+((tr3[,i])*1)
}


LUSC_C_te3$PI<-PI_te3
head(LUSC_C_te3)

#save selected data with PI value
#write.table(cbind("ID"=rownames(LUSC_C_te3), LUSC_C_te3),file="Ext_test_LTS_STS_with_11_genes_PI.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(LUSC_C_tr3), LUSC_C_tr3),file="Ext_test_LTS_MTS_STS_with_8genes_PI.txt",sep="\t",quote=F, row.names=F)


te_PI3 <- as.data.frame(LUSC_C_te3$PI)
rownames(te_PI3) <- rownames(LUSC_C_te3)
colnames(te_PI3) <- c("PI")
head(te_PI3)
##write.table(cbind("ID"=rownames(tr_PI3 ), tr_PI3 ),file="test_PI_based_on_11genes.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(tr_PI3 ), tr_PI3 ),file="test_PI_based_on_8genes.txt",sep="\t",quote=F, row.names=F)

######################################## Survival Object ############################################
surv_object_te3<- Surv(time = LUSC_C_te3$OS_month, event = LUSC_C_te3$Death_event)
surv_object_te3
dim(surv_object_te3 )


# ###survival analysis: fits cox ph model to find HR for PI
fit_te3 <- survfit(surv_object_te3~(te_PI3$PI>mean(te_PI3$PI)), data=te_PI3)
#fit_te <- survfit(surv_object_te~(LUSC_C_te$PI>mean(LUSC_C_te$PI)), data=LUSC_C_te)
summary(fit_te2)

#fit.coxph_te <- coxph(surv_object_te ~(LUSC_C_te$PI>mean(LUSC_C_te$PI)), data=LUSC_C_te)
fit.coxph_te3 <- coxph(surv_object_te3 ~(te_PI3$PI>mean(te_PI3$PI)), data=te_PI3)
summary(fit.coxph_te3)

fit.coxph_te$linear.predictors


te_res3 <- summary(fit.coxph_te3)
coeff3 <- round(te_res3$coefficients[1],2)
HR3 <- round(te_res3$conf.int[1],2)
int1_3 <- round(te_res3$conf.int[3],2)
int2_3 <- round(te_res3$conf.int[4],2)
CI3 <- round (te_res3$concordance[1],2)
pval3 <- te_res3$sctest[3]
pval1_3 <- format(pval3, scientific=T)
pval1_3

HR1_3 <- paste0(HR3 , " (",  int1_3,  " - ", int2_3, ") " )
HR1_3

te_res1 <- cbind(coeff3, HR1_3, CI3, pval1_3)
te_res1

names(te_res1) <- c("set","Beta_Coeff","HR with 95% CI","C-Index","P-value")
row.names(te_res1) <- c("Test_data")
#save results as a file
#write.table(cbind("ID"=rownames(tr_res1_3), tr_res1_3),file="Ext_test_results_based_on_11genes_PI.txt",sep="\t",quote=F, row.names=F)
#write.table(cbind("ID"=rownames(tr_res1_3), tr_res1_3),file="Ext_test_results_based_on_8genes_PI.txt",sep="\t",quote=F, row.names=F)

pp3 <-  ggsurvplot(fit_te3, data=te_PI3,
                   pval=TRUE,
                   risk.table=TRUE,
                   tables.height = 0.3, #add risk table & height
                   xlab="Time in Months",
                   risk.table.col="strata", break.time.by = 12,
                   conf.int = F, censor = TRUE,
                   #title= paste0("KM plot based on PI/Risk score"),
                   surv.median.line = "hv", # Add medians survival
                   palette=c("dodgerblue2", "red"), #add desired color
                   size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                   #legend.title = paste0(pathway),
                   legend.labs = c("Less than mean " , " more than Mean"),
                   font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
pp3


#save KM survival plot
#jpeg(file="KM_plot.jpeg", units="in", width=10, height=10, res=300)
pp3 <-  ggsurvplot(fit_te3, data=te_PI3,
                   #pval=TRUE,
                   risk.table=TRUE,
                   tables.height = 0.3, #add risk table & height
                   xlab="Time in Months",
                   risk.table.col="strata", break.time.by = 12,
                   conf.int = F, censor = TRUE,
                   #title= paste0("KM plot based on PI/Risk score"),
                   surv.median.line = "hv", # Add medians survival
                   palette=c("dodgerblue2", "red"), #add desired color
                   size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                   #legend.title = paste0(pathway),
                   legend.labs = c("Less than mean " , " more than Mean"),
                   font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
pp3

# customised the plot
pp3$plot <- pp3$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 1, hjust = 1,
    #label = "HR = 0.9 \n p < 0.001",
    #label =  paste0("HR = ", HR1, \n, "p-val <", pval),
    label =  paste0("HR = ", HR1_3, "\n",  "p-val = ", pval1_3,  "\n",  "C-Index = ", CI3),
    size = 5)

# now plot
pp3

jpeg("Ext_KM_plot_11genes_with_coeff.jpg", units="in", width=10, height=10, res=300)

#jpeg("Ext_KM_plot_8_genes.jpg", units="in", width=10, height=10, res=350)
print(pp3, newpage = FALSE)
dev.off()


#with quantiles
fit_te3_Q <- survfit(surv_object_te3~ cut(te_PI3$PI, quantile(te_PI3$PI)), data=te_PI3)
summary(fit_te3_Q)


fit.coxph_te3_Q <- coxph(surv_object_te3 ~ cut(te_PI3$PI, quantile(te_PI3$PI)), data=te_PI3)
summary(fit.coxph_te3_Q)




te_res3_Q <- summary(fit.coxph_te3_Q)
coeff3_Q <- round(te_res3_Q$coefficients[1],2)
HR3_Q <- round(te_res3_Q$conf.int[1],2)
int1_3_Q <- round(te_res3_Q$conf.int[3],2)
int2_3_Q <- round(te_res3_Q$conf.int[4],2)
CI3_Q <- round (te_res3_Q$concordance[1],2)
pval3_Q <- te_res3_Q$sctest[3]
pval1_3_Q <- format(pval3_Q, scientific=T)
pval1_3_Q

HR1_3_Q <- paste0(HR3_Q , " (",  int2_3_Q,  " - ", int1_3_Q, ") " )
HR1_3_Q

te_res1_3_Q <- cbind(coeff3_Q, HR1_3_Q, CI3_Q, pval1_3_Q)
te_res1_3_Q

######### KM plots with Quantiles ####

ggsurvplot(fit_te3_Q, data=te_PI3,
           #pval=teUE,
           risk.table=TRUE,
           xlab="Time in Months",
           tables.height = 0.35)

pp3_Q <- ggsurvplot(fit_te3_Q, data=te_PI3,
                    #pval=TRUE,
                    risk.table=TRUE,
                    xlab="Time in Months",
                    tables.height = 0.35,
                    surv.median.line = "hv", # Add medians survival
                    palette=c("green","dodgerblue2", "pink","red"), #add desired color
                    size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                    #title= paste0("KM plot based on Quantiles of PI/Risk score"),
                    legend.title = paste0("Quantiles of PI score"),
                    #legend.title = paste0(pathway),
                    legend.labs = c("Q1" , "Q2", "Q3", "Q4"),
                    risk.table.col="strata", break.time.by = 6)
#fit_tr <- survfit(surv_object_tr~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
pp3_Q 

# customised the plot
pp3_Q$plot <- pp3_Q$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 1, hjust = 1,
    #label = "HR = 0.9 \n p < 0.001",
    #label =  paste0("HR = ", HR1, \n, "p-val <", pval),
    label =  paste0("HR = ", HR1_3_Q, "\n",  "p-val = ", round(pval3_Q,3),  "\n",  "C-Index = ", CI3_Q),
    size = 5)

# now plot
pp3_Q



#jpeg("Ext_KM_plot_11genes_with_coeff_with_Quantiles.jpg", units="in", width=10, height=10, res=300)
jpeg("Ext_KM_plot_8genes_with_coeff_with_Quantiles.jpg", units="in", width=10, height=10, res=300)

print(pp3_Q, newpage = FALSE)
dev.off()

####### ROC plots 

##Create ROC plot for 1-,3- and 5-years survival time prediction
te_roc1 <- survivalROC(Stime        = LUSC_C_te3$OS_month,
                       status       = LUSC_C_te3$Death_event,
                       marker       = LUSC_C_te3$PI,
                       predict.time = 12,
                       method       = "KM", 
                       #lambda = lambda_min , 
                       span = NULL,
                       window ="symmetric")
te_roc1

te_roc3 <- survivalROC(Stime        = LUSC_C_te3$OS_month,
                       status       = LUSC_C_te3$Death_event,
                       marker       = LUSC_C_te3$PI,
                       predict.time = 36,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmeteic")
te_roc3

te_roc4 <- survivalROC(Stime        = LUSC_C_te3$OS_month,
                       status       = LUSC_C_te3$Death_event,
                       marker       = LUSC_C_te3$PI,
                       predict.time = 48,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
te_roc4
te_roc5 <- survivalROC(Stime        = LUSC_C_te3$OS_month,
                       status       = LUSC_C_te3$Death_event,
                       marker       = LUSC_C_te3$PI,
                       predict.time = 60,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
te_roc5
te_roc6 <- survivalROC(Stime        = LUSC_C_te3$OS_month,
                       status       = LUSC_C_te3$Death_event,
                       marker       = LUSC_C_te3$PI,
                       predict.time = 72,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")



plot(te_roc1$FP, te_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red",
     ylab="TP",main= "AUC Curve for Survival Prediction")
lines(te_roc3$FP, te_roc3$TP, type="l", lty=2, col="blue")
lines(te_roc5$FP, te_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(te_roc1$AUC,2)),  paste("3 Years AUC = ",round(te_roc3$AUC,2)),  paste("5 Years AUC = ",round(te_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)


jpeg(file="ROC_test_11genes_with_coeff.jpeg", units="in", width=10, height=10, res=300)

#jpeg(file="ROC_test_8genes_with_coeff.jpeg", units="in", width=10, height=10, res=300)
plot(te_roc1$FP, te_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red",
     ylab="TP",main= "AUC Curve for Survival Prediction")
lines(te_roc3$FP, te_roc3$TP, type="l", lty=2, col="blue")
lines(te_roc5$FP, te_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(te_roc1$AUC,2)),  paste("3 Years AUC = ",round(te_roc3$AUC,2)),  paste("5 Years AUC = ",round(te_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)
dev.off()




######### multivariate #########

data1_new <- read.table("train_8_genes_with_PI_with_Clin_data.txt", sep="\t", header=T, row.names=1, check.names = F)
# create survival object
surv_object_tr <- Surv(time = data1_new$OS_month, event = data1_new$Death_Status)#Create survival object


cox_multivariate_tr1 <- coxph(surv_object_tr ~  Gender + as.numeric(Age) + IDH_mutation_status +  
                                MGMTp_meth_status + Codel_1p19q_status ,  data=data1_new )

cox_multivariate_tr1 <- coxph(surv_object_tr ~  Gender + as.numeric(Age) + IDH_mutation_status +  
                                MGMTp_meth_status  ,  data=data1_new )

summary(cox_multivariate_tr1 )
ggforest(cox_multivariate_tr1, data=data1_new )

#11 genes
cox_multivariate_tr2 <- coxph(surv_object_tr ~  Gender + as.numeric(Age) + IDH_mutation_status +  
                                MGMTp_meth_status +
                                CPZ + FBXW2 + FCGR3A + KBTBD2 + MYO1B + PMM2 + PRKCI + 
                                SSNA1 + TINF2 + TNC + UTP20,  data=data1_new )

#8 genes
cox_multivariate_tr2 <- coxph(surv_object_tr ~  Gender + as.numeric(Age) + IDH_mutation_status +  
                                MGMTp_meth_status +
                                CPZ + FBXW2 + PMM2 + PRKCI +  MYO1B +
                                SSNA1 + TNC + UTP20,  data=data1_new )


summary(cox_multivariate_tr2 )
ggforest(cox_multivariate_tr2, data=data1_new )

jpeg(file="multivariate_train_with_11genes.jpeg", units="in", width=10, height=10, res=350)
jpeg(file="multivariate_train_with_8genes.jpeg", units="in", width=10, height=10, res=350)
ggforest(cox_multivariate_tr2, data=data1_new, fontsize = 1.0 )
dev.off()


cox_multivariate_tr3 <- coxph(surv_object_tr ~  Gender + as.numeric(Age) + IDH_mutation_status +  
                                MGMTp_meth_status + Codel_1p19q_status +PI ,  data=data1_new)

cox_multivariate_tr3 <- coxph(surv_object_tr ~  Gender + as.numeric(Age) + IDH_mutation_status +  
                                MGMTp_meth_status +PI ,  data=data1_new)

summary(cox_multivariate_tr3 )
ggforest(cox_multivariate_tr3, data=data1_new )

jpeg(file="multivariate_train_with_PI.jpeg", units="in", width=10, height=10, res=350)

jpeg(file="multivariate_train_with_8genes_PI.jpeg", units="in", width=10, height=10, res=350)
ggforest(cox_multivariate_tr3, data=data1_new, fontsize = 1.0 )
dev.off()

data_test1_new <- read.table("Ext_test_LTS_STS_with_11_genes_PI_with_Clin_data.txt", sep="\t", header=T, row.names=1, check.names = F)
data_test1_new <- read.table("Ext_test_LTS_STS_with_8_genes_PI_with_Clin_data.txt", sep="\t", header=T, row.names=1, check.names = F)

data_test1_new
surv_object_te <- Surv(time = data_test1_new$OS_month, event = data_test1_new$Death_event)#Create survival object

head(data_test1_new [1:15])

cox_multivariate_te1 <- coxph(surv_object_te ~  Gender + as.numeric(Age) + IDH_mutation_status +  
                                MGMTp_meth_status + Codel_1p19q_status + PI ,  data=data_test1_new)
cox_multivariate_te1 <- coxph(surv_object_te ~  Gender + as.numeric(Age) + IDH_mutation_status +  
                                MGMTp_meth_status  ,  data=data_test1_new)


summary(cox_multivariate_te1 )
ggforest(cox_multivariate_te1, data=data_test1_new, fontsize = 1.0)


# 11 genes
cox_multivariate_te2 <- coxph(surv_object_te ~  Gender + as.numeric(Age) + IDH_mutation_status +  
                                MGMTp_meth_status + 
                                CPZ + FBXW2 + FCGR3A + KBTBD2 + MYO1B + PMM2 + PRKCI + 
                                SSNA1 + TINF2 + TNC + UTP20,  data=data_test1_new )

# 8 genes
cox_multivariate_te2 <- coxph(surv_object_te ~  Gender + as.numeric(Age) + IDH_mutation_status +  
                                MGMTp_meth_status + 
                                CPZ + FBXW2 +  MYO1B + PMM2 + PRKCI + 
                                SSNA1 +  TNC + UTP20,  data=data_test1_new )



summary(cox_multivariate_te2 )
ggforest(cox_multivariate_te2, data=data_test1_new , fontsize = 1.0)

jpeg(file="multivariate_test_with_11_genes.jpeg", units="in", width=10, height=10, res=350)

jpeg(file="multivariate_test_with_8_genes.jpeg", units="in", width=10, height=10, res=350)
ggforest(cox_multivariate_te2, data=data_test1_new, fontsize = 1.0 )
dev.off()


cox_multivariate_te3 <- coxph(surv_object_te ~  Gender + as.numeric(Age) + IDH_mutation_status +  
                                MGMTp_meth_status + PI ,  data=data_test1_new)


summary(cox_multivariate_te3 )
ggforest(cox_multivariate_te3, data=data_test1_new )

jpeg(file="multivariate_test_with_PI.jpeg", units="in", width=10, height=10, res=350)

jpeg(file="multivariate_test_with_8genes_PI.jpeg", units="in", width=10, height=10, res=350)
ggforest(cox_multivariate_te3, data=data_test1_new, fontsize = 1.0 )
dev.off()



########################################  Nomogram   ########################################## 

# with 11 genes #######
#tr_clin_new <- read.table("tr_PI_data_with_clin.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
tr_clin_new <- read.table("train_11_genes_with_PI_with_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
te_clin_new <- read.table("Ext_test_LTS_STS_with_11_genes_PI_with_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)


# with 8 genes 
tr_clin_new <- read.table("train_8_genes_with_PI_with_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
te_clin_new <- read.table("Ext_test_LTS_STS_with_8_genes_PI_with_Clin_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)


tr_clin_new$Codel_1p19q_status

head(tr_clin_new,2)
dim(tr_clin_new)

# create survival object
surv_object_p_tr <- Surv(time = as.numeric(tr_clin_new$OS_month), event = tr_clin_new$Death_Status)
surv_object_te <- Surv(time = te_clin_new$OS_month, event = te_clin_new$Death_event)#Create survival object


d <-  cbind(tr_clin_new["OS_month"], tr_clin_new["Death_Status"], tr_clin_new["Gender"],tr_clin_new["Age"], tr_clin_new["IDH_mutation_status"],  tr_clin_new["MGMTp_meth_status"],  tr_clin_new["Codel_1p19q_status"],tr_clin_new["PI"])
d_te <-  cbind(te_clin_new["OS_month"], te_clin_new["Death_event"],  te_clin_new["Gender"],te_clin_new["Age"], te_clin_new["IDH_mutation_status"], te_clin_new["MGMTp_meth_status"], te_clin_new["Codel_1p19q_status"], te_clin_new["PI"])

dim(d)
ddist <- datadist(d)
options(datadist='ddist')


d
ddist_te <- datadist(d_te)
options(datadist='ddist_te')

### logistic regression models
f2 <- lrm(Death_Status ~ Age + Gender +  IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status + PI , x = T,y = T, data = d)
f <- lrm(Death_Status ~ Age + Gender +  IDH_mutation_status +  MGMTp_meth_status + PI , x = T,y = T, data = d)

cal_f <- calibrate(f, method ="boot", B = 50) 
plot(cal_f)
plot(cal_f,xlim = c(0,1.0),ylim = c(0,1.0))

f3 <- lrm(Death_Status ~ Age +   IDH_mutation_status +  PI , x = T,y = T, data = d)

cal_f3 <- calibrate(f3, method ="boot", B = 50) 
plot(cal_f3)
plot(cal_f3,xlim = c(0,1.0),ylim = c(0,1.0))

f1 <- lrm(Death_Status ~  PI , data = d,  x = T,y = T)



#nom <- nomogram(f, fun=plogis, funlabel="Risk of Death")
nom <- nomogram(f, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
nom1 <- nomogram(f1, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")
nom2 <- nomogram(f2, fun= function(x)1/(1+exp(-x)),  lp = F, funlabel = "Risk")


plot(nom)

jpeg("Nomogram_train_Risk_LR.jpg", units="in", width=10, height=10, res=300)
plot(nom)
dev.off()


jpeg("Nomogram_PI_Grade_IDH_codel_Risk_LR.jpg", units="in", width=10, height=10, res=300)
plot(nom4)
dev.off()

jpeg("Nomogram_PI_Gender_IDH_MGMT_age_Risk_LR.jpg", units="in", width=10, height=10, res=300)
plot(nom5)
dev.off()

f
Survival(f)
surv1<- Survival(f)

#coxph models
cox1 <-cph(Surv(OS_month,Death_Status==1) ~ Age + Gender +  IDH_mutation_status + MGMTp_meth_status + Codel_1p19q_status  + PI,  x = T,y = T, data = d, surv = T)


##### #################
surv<- Survival(cox1)
risk <- function(x)1/(1+exp(-x))
surv_1<- function(x)surv(1*12,lp=x) # defined time.inc,1 year OS
surv_2<- function(x)surv(1*36,lp=x) # defined time.inc,3 year OS
surv_3<- function(x)surv(1*60,lp=x) # defined time.inc,5 year OS

nom_cox1<-nomogram(cox1,fun = list(risk, surv_1,surv_2,surv_3),lp = F,
                   funlabel = c("Risk", "1-Year Survival Probability", "3-Year Survival Probability","5-Year Survival Probability"),
                   maxscale = 100,
                   fun.at = c('1.0','0.95','0.90','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom_cox1),xfrac = .3)

jpeg("Nomogram_train_PI_OS_time_COX.jpg", units="in", width=15, height=10, res=350)
jpeg("Nomogram_train_8genes_based_PI_OS_time_COX1.jpg", units="in", width=15, height=10, res=350)
plot(nom_cox1, xfrac = .3)
dev.off()


v <- validate(cox1 , dxy = TRUE, B = 1000) 
Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"] 
orig_Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.orig"] 
bias_corrected_c_index  <- abs(Dxy)/2+0.5 
orig_c_index <- abs(orig_Dxy)/2+0.5  
bias_corrected_c_index
orig_c_index



###########

cox <-cph(Surv(OS_month,Death_Status==1) ~ Age + Gender +  IDH_mutation_status + MGMTp_meth_status  + PI,  x = T,y = T, data = d, surv = T)

v <- validate(cox , dxy = TRUE, B = 1000) 
Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"] 
orig_Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.orig"] 
bias_corrected_c_index  <- abs(Dxy)/2+0.5 
orig_c_index <- abs(orig_Dxy)/2+0.5  
bias_corrected_c_index
orig_c_index


##### #################
surv<- Survival(cox)
risk <- function(x)1/(1+exp(-x))
surv_1<- function(x)surv(1*12,lp=x) # defined time.inc,1 year OS
surv_2<- function(x)surv(1*36,lp=x) # defined time.inc,3 year OS
surv_3<- function(x)surv(1*60,lp=x) # defined time.inc,5 year OS

nom_cox<-nomogram(cox,fun = list(risk, surv_1,surv_2,surv_3),lp = F,
                  funlabel = c("Risk", "1-Year Survival Probability", "3-Year Survival Probability","5-Year Survival Probability"),
                  maxscale = 100,
                  fun.at = c('1.0','0.95','0.90','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom_cox),xfrac = .3)

jpeg("Nomogram_train_PI_OS_time_COX.jpg", units="in", width=15, height=10, res=350)
jpeg("Nomogram_train_8genes_based_PI_OS_time_COX.jpg", units="in", width=15, height=10, res=350)
plot(nom_cox, xfrac = .3)
dev.off()

###########################
cox_1 <-cph(Surv(OS_month,Death_Status==1) ~ PI,  x = T,y = T, data = d, surv = T)

cox2 <-cph(Surv(OS_month,Death_Status==1) ~ Age +  Gender +  IDH_mutation_status + Codel_1p19q_status + PI,  x = T,y = T, data = d, surv = T)


v <- validate(cox2 , dxy = TRUE, B = 1000) 
Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"] 
orig_Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.orig"] 
bias_corrected_c_index  <- abs(Dxy)/2+0.5 
orig_c_index <- abs(orig_Dxy)/2+0.5  
bias_corrected_c_index
orig_c_index


surv2<- Survival(cox2)
risk2 <- function(x)1/(1+exp(-x))
surv_1_2<- function(x)surv2(1*12,lp=x) # defined time.inc,1 year OS
surv_2_2<- function(x)surv2(1*36,lp=x) # defined time.inc,3 year OS
surv_3_2<- function(x)surv2(1*60,lp=x) # defined time.inc,5 year OS

nom_cox2<-nomogram(cox2,fun = list(risk2, surv_1_2,surv_2_2,surv_3_2),lp = F,
                   funlabel = c("Risk", "1-Year Survival Probability", "3-Year Survival Probability","5-Year Survival Probability"),
                   maxscale = 100,
                   fun.at = c('1.0','0.95','0.90','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom_cox2),xfrac = .9)

jpeg("Nomogram_train_8_genes_based_PI_OS_time_COX2.jpg", units="in", width=15, height=10, res=350)
plot(nom_cox2, xfrac = .3)
dev.off()



###################
cox3 <-cph(Surv(OS_month,Death_Status==1) ~ Age +  IDH_mutation_status +   PI,  x = T,y = T, data = d, surv = T)


v <- validate(cox3 , dxy = TRUE, B = 1000) 
Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.corrected"] 
orig_Dxy = v[rownames(v)=="Dxy", colnames(v)=="index.orig"] 
bias_corrected_c_index  <- abs(Dxy)/2+0.5 
orig_c_index <- abs(orig_Dxy)/2+0.5  
bias_corrected_c_index
orig_c_index

#################
surv3<- Survival(cox3)
risk3 <- function(x)1/(1+exp(-x))
surv_1_3<- function(x)surv3(1*12,lp=x) # defined time.inc,1 year OS
surv_2_3<- function(x)surv3(1*36,lp=x) # defined time.inc,3 year OS
surv_3_3<- function(x)surv3(1*60,lp=x) # defined time.inc,5 year OS

nom_cox3<-nomogram(cox3,fun = list(risk3, surv_1_3,surv_2_3,surv_3_3),lp = F,
                   funlabel = c("Risk", "1-Year Survival Probability", "3-Year Survival Probability","5-Year Survival Probability"),
                   maxscale = 100,
                   fun.at = c('1.0','0.95','0.90','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom_cox3),xfrac = .9)

jpeg("Nomogram_train_PI_OS_time_COX3.jpg", units="in", width=15, height=10, res=350)
jpeg("Nomogram_train_8_genes_based_PI_OS_time_COX3.jpg", units="in", width=15, height=10, res=350)
plot(nom_cox3, xfrac = .3)
dev.off()




###### calibration curve ####
cal1<- calibrate(cox3, cmethod = "KM", method = "boot", u = 12, m = 20, B = 100)
jpeg("Calibration_plot_cox3_1yrs .jpg", units="in", width=15, height=10, res=350)
plot(cal1,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 1-Year",
     ylab="Actual Probability of 1-Year",
     col=c(rgb(192,98,83,maxColorValue=255)))
dev.off()

cal2<- calibrate(cox3, cmethod = "KM", method = "boot", u = 24, m = 20, B = 100)
jpeg("Calibration_plot_cox3_2yrs .jpg", units="in", width=15, height=10, res=350)
plot(cal2,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 2-Year",
     ylab="Actual Probability of 2-Year",
     col=c(rgb(192,98,83,maxColorValue=255)))
dev.off()


cal3<- calibrate(cox3, cmethod = "KM", method = "boot", u = 36, m = 25, B = 100)
jpeg("Calibration_plot_cox3_3yrs .jpg", units="in", width=15, height=10, res=350)
plot(cal3,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 3-Year",
     ylab="Actual Probability of 3-Year",
     col=c(rgb(192,98,83,maxColorValue=255)))
dev.off()
cal4<- calibrate(cox3, cmethod = "KM", method = "boot", u = 48, m = 30, B = 100)
plot(cal4,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 4-Year",
     ylab="Actual Probability of 4-Year",
     col=c(rgb(192,98,83,maxColorValue=255)))



cal5<- calibrate(cox3, cmethod = "KM", method = "boot", u = 60, m = 20, B = 100)

jpeg("Calibration_plot_cox3_5yrs .jpg", units="in", width=15, height=10, res=350)
plot(cal5,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 5-Year",
     ylab="Actual Probability of 5-Year",
     col=c(rgb(192,98,83,maxColorValue=255)))
dev.off()


##################
###### calibration curve ####
cal1<- calibrate(cox2, cmethod = "KM", method = "boot", u = 12, m = 20, B = 100)
jpeg("Calibration_plot_cox3_1yrs .jpg", units="in", width=15, height=10, res=350)
jpeg("Calibration_plot_8genes_cox2_1yrs .jpg", units="in", width=15, height=10, res=350)
plot(cal1,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 1-Year",
     ylab="Actual Probability of 1-Year",
     col=c(rgb(192,98,83,maxColorValue=255)))
dev.off()

cal2<- calibrate(cox2, cmethod = "KM", method = "boot", u = 24, m = 20, B = 100)
jpeg("Calibration_plot_8_genes_PI_cox2_2yrs .jpg", units="in", width=15, height=10, res=350)
plot(cal2,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 2-Year",
     ylab="Actual Probability of 2-Year",
     col=c(rgb(192,98,83,maxColorValue=255)))
dev.off()


cal3<- calibrate(cox2, cmethod = "KM", method = "boot", u = 36, m = 20, B = 100)
jpeg("Calibration_plot_8genes_cox2_3yrs .jpg", units="in", width=15, height=10, res=350)
plot(cal3,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 3-Year",
     ylab="Actual Probability of 3-Year",
     col=c(rgb(192,98,83,maxColorValue=255)))
dev.off()
cal4<- calibrate(cox2, cmethod = "KM", method = "boot", u = 48, m = 20, B = 100)
plot(cal4,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 4-Year",
     ylab="Actual Probability of 4-Year",
     col=c(rgb(192,98,83,maxColorValue=255)))



cal5<- calibrate(cox2, cmethod = "KM", method = "boot", u = 60, m = 20, B = 100)

jpeg("Calibration_plot_8genes_cox2_5yrs .jpg", units="in", width=15, height=10, res=350)
plot(cal5,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 5-Year",
     ylab="Actual Probability of 5-Year",
     col=c(rgb(192,98,83,maxColorValue=255)))
dev.off()



#################
sum.surv_tr2<-summary(cox)

sum.surv_tr2

c_index_tr2<-sum.surv_tr2$
  
  c_index_tr2

sum.surv_tr2$sctest

###################
surv1<- Survival(cox1)
risk_1 <- function(x)1/(1+exp(-x))
surv1_1<- function(x)surv1(1*12,lp=x) # defined time.inc,1 year OS
surv1_2<- function(x)surv1(1*36,lp=x) # defined time.inc,3 year OS
surv1_3<- function(x)surv1(1*60,lp=x) # defined time.inc,5 year OS

nom1_cox<-nomogram(cox1,fun = list(risk_1, surv1_1,surv1_2,surv1_3),lp = F,
                   funlabel = c("Risk", "1-Year Survival Probability", "3-Year Survival Probability","5-Year Survival Probability"),
                   maxscale = 100,
                   fun.at = c('0.95','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom1_cox),xfrac = .7)

jpeg("Nomogram_PI_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
plot(nom1_cox)
dev.off()


########## C-index and p-value ######

#cox1$coefficients
#f<-coxph(Surv(OS_month,Death_Status==1) ~  Age + Gender +  IDH_mutation_status  + MGMTp_meth_status + PI , data = d)
cox_te <-cph(Surv(OS_month,Death_event==1) ~ Age + Gender +  IDH_mutation_status + MGMTp_meth_status  + PI,  x = T,y = T, data = d_te, surv = T)
cox_te
###################
surv1_te <- Survival(cox_te)
risk_1_te <- function(x)1/(1+exp(-x))
surv1_1_te<- function(x)surv1_te(1*12,lp=x) # defined time.inc,1 year OS
surv1_2_te<- function(x)surv1_te(1*36,lp=x) # defined time.inc,3 year OS
surv1_3_te<- function(x)surv1_te(1*60,lp=x) # defined time.inc,5 year OS

nom1_cox_te<-nomogram(cox_te,fun = list(risk_1_te, surv1_1_te,surv1_2_te,surv1_3_te),lp = F,
                      funlabel = c("Risk", "1-Year Survival Probability", "3-Year Survival Probability","5-Year Survival Probability"),
                      maxscale = 100,
                      fun.at = c('0.95', '0.9','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1')  )
plot((nom1_cox_te),xfrac = .9)

jpeg("Nomogram_test_PI_OS_time_COX.jpg", units="in", width=15, height=10, res=300)
plot((nom1_cox_te),xfrac = .9)
dev.off()


########## concordance Index ######
res_con <- rcorrcens(Surv(OS_month,Death_Status) ~ predict(cox3), data = tr_clin_new)
res_con[1]

Cindex = 1- res_con[1]
Cindex
###############################

################### Test data ########
########## C-index and p-value ######


#Prediction on Test data
#f_c_te<-coxph(Surv(OS_month,Death_event==1) ~ PI + Age + Gender +  IDH_mutation_status + MGMTp_meth_status , data = d_te)

Surv.obj_test=with(d_te,Surv(OS_month,Death_event))
###Create your survival estimates

estimates_1=survest(cox1,newdata=d_te,times=12)$surv
estimates_3=survest(cox1,newdata=d_te,times=36)$surv
estimates_5=survest(cox1,newdata=d_te,times=60)$surv

estimates_1=survest(cox3,newdata=d_te,times=12)$surv
estimates_3=survest(cox3,newdata=d_te,times=36)$surv
estimates_5=survest(cox3,newdata=d_te,times=60)$surv


estimates_1=survest(cox2,newdata=d_te,times=12)$surv
estimates_3=survest(cox2,newdata=d_te,times=36)$surv
estimates_5=survest(cox2,newdata=d_te,times=60)$surv
#estimates

###Determine concordance
rcorr.cens(x=estimates_1,S=Surv.obj_test)
rcorr.cens(x=estimates_3,S=Surv.obj_test)
rcorr.cens(x=estimates_5,S=Surv.obj_test)

############### Prediction of samples ########
#Predict(f, Age=43,  Gender='Male', IDH_mutation_status="Mutant", MGMTp_meth_status="methylated", PI=3.3, fun=plogis)

Predict(cox, Age=43,  Gender='Male', IDH_mutation_status="Mutant", MGMTp_meth_status="methylated", PI=3.3, fun=plogis)
Predict(cox3, Age=53,   IDH_mutation_status="Wildtype", PI=3.35, fun=plogis)
Predict(cox3, Age=42,   IDH_mutation_status="Wildtype", PI=3.76, fun=plogis)
Predict(cox3, Age=53,   IDH_mutation_status="Wildtype", PI=3.35, fun=plogis)
Predict(cox3, Age=79,   IDH_mutation_status="Wildtype", PI=4.01, fun=plogis)
Predict(cox2, Age=24,   IDH_mutation_status="Mutant", PI=3.477, fun=plogis)

fit3 <- lrm(Death_Status ~  Age + Gender +  PI , data = d)

d$MGMTp_meth_status
# type of predicted value
predict(fit3,type="lp") #linear predictor ("lp")
predict(fit3,type="expected") #expected number of events given the covariates and follow-up time 
predict(fit3,type="risk",se.fit=TRUE) # risk score exp(lp)
predict(fit3,type="terms",se.fit=TRUE) #terms of the linear predictor
predict(fit3, newdata=te_clin_new , type="survival",se.fit=TRUE) #terms of the linear predictor

######################### ROC plots  #######

library(ROCR) 
cox7$sformula
cox2$linear.predictors
predvalue2 <- predict(cox)
pred2 <- prediction(predvalue2 , d$Death_Status)
pred2
auc <- performance(pred2,"auc")
auc

perf<- performance(pred2,"tpr","fpr")
perf
plot(perf)
abline(0,1, col = 3, lty = 2)
##Create ROC plot for 1-,3- and 5-years survival time

tr_roc1 <- survivalROC(Stime        = d$OS_month,
                       status       = d$Death_Status,
                       marker       = predvalue2,
                       predict.time = 12,
                       method       = "KM", 
                       #lambda = lambda_min , 
                       span = NULL,
                       window ="symmetric")
tr_roc1

tr_roc3 <- survivalROC(Stime        = d$OS_month,
                       status       = d$Death_Status,
                       marker       = predvalue2,
                       predict.time = 48,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc3
tr_roc5 <- survivalROC(Stime        = d$OS_month,
                       status       = d$Death_Status,
                       marker       = predvalue2,
                       predict.time = 60,
                       method       = "KM", 
                       #lambda = lambda_min, 
                       span = NULL, 
                       window ="symmetric")
tr_roc5


tr_roc6
jpeg(file="ROC_train_PI_Grade.jpeg", units="in", width=10, height=10, res=300)
plot(tr_roc1$FP, tr_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red",
     ylab="TP",main= "AUC Curve for Survival Prediction based on PI+Grade Model")
lines(tr_roc3$FP, tr_roc3$TP, type="l", lty=2, col="blue")
lines(tr_roc5$FP, tr_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("4 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)
dev.off()




########################## Calibration curves #############

?calibrate



cal1 <- calibrate(cox, cmethod="KM", method="boot", u = 20,  B = 100)
plot(cal1)
cal1 <- calibrate(cox, cmethod="KM", method="boot", u = 20,  B = 100)
plot(cal1)
cal1 <- calibrate(cox, cmethod="KM", method="boot",  u = 365, m = 50, B = 100)
plot(cal1)
plot(cal1,lwd=2,lty=1,
     errbar.col=c(rgb(0,118,192,maxColorValue=255)),
     xlim=c(0.1,1),ylim=c(0.1,1),
     xlab="Nomogram-Predicted Probability of 1-Year",
     ylab="Actual Probability of 1-Year",
     col=c(rgb(192,98,83,maxColorValue=255)))


surv1<- Survival(cox1)
surv2<- Survival(cox2)
#surv1_1<- function(x)surv(1*12,lp=x) # defined time.inc,1 year OS
#surv2_2<- function(x)surv(1*36,lp=x) # defined time.inc,3 year OS
#surv3_3<- function(x)surv(1*60,lp=x) # defined time.inc,5 year OS




lrm(Death_Status ~  Grade +  IDH_mutation_status + Codel_1p19q_status + PI , data = d)

head(d)
f2 <- psm(Surv(OS_month,Death_Status) ~  Grade +  IDH_mutation_status + Codel_1p19q_status +PI, data=d, dist='lognormal')
med  <- Quantile(f2)
surv1 <- Survival(f2)  # This would also work if f was from cph
surv1
plot(nomogram(f2, fun=function(x) med(lp=x), funlabel="Median Survival Time"))
nom2 <- nomogram(f2, fun=list(function(x) surv1(12, x),
                              function(x) surv1(36, x),
                              function(x) surv1(60, x),
                              function(x) surv1(120, x)),
                 funlabel=c("12-Month Survival Probability", 
                            "36-month Survival Probability",
                            "60-month Survival Probability",
                            "120-month Survival Probability"))
plot(nom2, xfrac=.7)

f2$score
f2$dist
f2$non.slopes
surv1(12, x)

f2$scale.pred

cox2$score


###########################
surv_object3 <- Surv(time = data3$OS_month, event = data3$Death_Status)
fit3 <- coxph(surv_object3 ~ Grade  +PI, data=data3)

fit3 <- lrm(Death_Status ~  Age + Gender +  IDH_mutation_status + MGMTp_meth_status + PI , data = d)
fit3 <- lrm(Death_Status ~  Age + Gender +  PI , data = d)

d$MGMTp_meth_status
# type of predicted value
predict(fit3,type="lp") #linear predictor ("lp")
predict(fit3,type="expected") #expected number of events given the covariates and follow-up time 
predict(fit3,type="risk",se.fit=TRUE) # risk score exp(lp)
predict(fit3,type="terms",se.fit=TRUE) #terms of the linear predictor
predict(fit3, newdata=new_te, type="survival",se.fit=TRUE) #terms of the linear predictor



tr_CI<-concordance.index(x=train$Ridge.Classifier, surv.time=train$DFS.time, surv.event=train$DFS, method="noether")
#### x=predictions

tr_CI

D_index_tr<-D.index(x=train$pred1, surv.time=train$DFS.time, surv.event=train$DFS, alpha = 0.05,
                    method.test = c("logrank"), na.rm = T)
D_index_tr

HR_tr<-hazard.ratio(x=train$Ridge.Classifier, surv.time=train$DFS.time, surv.event=train$DFS, alpha = 0.05,
                    method.test = c("logrank"), na.rm = T)

HR_tr



plot(survfit(f, newdata=d,  xscale=365.25, xlab="Years", ylab="Survival", conf.int=F) )
# also plot the predicted survival for a 70 year old
lines(survfit(fit3, newdata=new_data), xscale=365.25, xlab="Years", ylab="Survival") 

