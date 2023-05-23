#Load required libraries
library(caret)
library(ggfortify)
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(ggfortify)

set.seed(7)
setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Primary_grade4_survival/")



data1 <- read.table("Primary_grade4_clin_data_STS_MTS_LTS_both_trts_with_exp_data.txt", header=T, sep="\t", row.names=1, check.names=F)



######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
sel_features_results<-read.table("3_features.txt",header =TRUE, sep = "\t", row.names=1, check.names = FALSE)
head(sel_features_results)
dim(sel_features_results)


#data preparation with selected features 
#sel_train <- as.data.frame(data1[,colnames(data1) %in% c(row.names(features)), ])
#sel_train <- as.data.frame(data1[,colnames(data1) %in% c(sel_features_results$ID), ])
sel_train <- as.data.frame(data1[,colnames(data1) %in% c(row.names(sel_features_results)), ])
head(sel_train,2)
dim(sel_train)


######### Rank Data #####
#ex <- read.table("example.txt",  header=T, sep="\t", row.names=1, check.names=F)
ex <- read.table("3genes_train.txt",  header=T, sep="\t", row.names=1, check.names=F)

ex 
tmp_list <- list()
for(i in seq(from=1, to=length(ex), by=1))
{
  
  ex1 <- ex[i]
  
  col1_r <-  apply(ex1, 2, rank)
  col1_r 
  
  tmp_list[[i]] <- col1_r 
  
}

df<- data.frame(tmp_list)

df

write.table(df,file="3_genes_rank_data.txt", sep='\t',  quote = F,row.names = TRUE)

sel_train<- read.table("3_genes_rank_data_tpose", header=T, sep="\t", row.names=1, check.names=F)

head(sel_train)

######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
train_feature_mat<- cbind(data1["OS_month"],data1["Death_Status"],sel_train)
head(train_feature_mat,2)
dim(train_feature_mat)

############ remove where OS.time=NA ############
train_feature_mat1<-subset(train_feature_mat,OS_month!="NA")
train_feature_mat1<-subset(train_feature_mat1,OS_month!=0)


#train_feature_mat2<- cbind(data2["OS_month"],data2["Group"],sel_train)
# save files with selected genes & survival information #########
write.table(cbind("ID"=rownames(train_feature_mat1), train_feature_mat1),file="sel_train.txt",sep="\t",quote=F, row.names=F)

#Create prognostic index for pathway
LUSC_C_tr=train_feature_mat1
tr <- LUSC_C_tr[3:ncol(LUSC_C_tr)]

E= length(tr)
E

head(tr,4)
tr[,2]
sel_features_results[2,1]
head(sel_features_results,2)

PI_tr=0
for(i in seq(from=1, to=E ,by=1))
{
  PI_tr= PI_tr+((tr[,i])*(sel_features_results[i,1]))
}



head(tr,3)
head(sel_features_results)

LUSC_C_tr$PI<-PI_tr

data1$RSI <- LUSC_C_tr$PI

df <- cbind(data1$Survival_class, data1$RSI)
df
#Find the quartiles (25th, 50th, and 75th percentiles) of the vector
quantile(data1$RSI, probs = c(.25, .5, .75))
head(LUSC_C_tr)

#save selected data with PI value
write.table(cbind("ID"=rownames(LUSC_C_tr), LUSC_C_tr),file="train_with_PI.txt",sep="\t",quote=F, row.names=F)

tr_PI <- as.data.frame(LUSC_C_tr$PI)
rownames(tr_PI) <- rownames(LUSC_C_tr)
colnames(tr_PI) <- c("PI")
write.table(cbind("ID"=rownames(tr_PI ), tr_PI ),file="tr_PI",sep="\t",quote=F, row.names=F)


######################################## Survival Object ############################################
surv_object_tr <- Surv(time = LUSC_C_tr$OS_month, event = LUSC_C_tr$Death_Status)
surv_object_tr
dim(surv_object_tr )


dim(tr)
head(tr_PI)
mean(tr_PI$PI)

# ###survival analysis: fits cox ph model to find HR for PI
fit_tr <- survfit(surv_object_tr~(tr_PI$PI>mean(tr_PI$PI)), data=tr_PI)
#fit_tr <- survfit(surv_object_tr~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
summary(fit_tr)

#fit.coxph_tr <- coxph(surv_object_tr ~(LUSC_C_tr$PI>mean(LUSC_C_tr$PI)), data=LUSC_C_tr)
fit.coxph_tr <- coxph(surv_object_tr ~(tr_PI$PI>mean(tr_PI$PI)), data=tr_PI)
summary(fit.coxph_tr)

fit.coxph_tr$linear.predictors


tr_res <- summary(fit.coxph_tr)
coeff <- round(tr_res$coefficients[1],2)
HR <- round(tr_res$conf.int[1],2)
int1 <- round(tr_res$conf.int[3],2)
int2 <- round(tr_res$conf.int[4],2)
CI <- round (tr_res$concordance[1],2)
pval <- tr_res$sctest[3]

HR1 <- paste0(HR , " (",  int1,  " - ", int2, ") " )
HR1

tr_res1 <- cbind(coeff, HR1, CI, pval)
tr_res1
#save results as a file
write.table(cbind("ID"=rownames(tr_res1), tr_res1),file="tr_res1.txt",sep="\t",quote=F, row.names=F)


#save KM survival plot
#jpeg(file= paste0(path, subfolder[j], "/", "KM_plot_train_with_PI.jpeg"), units="in", width=10, height=10, res=300)
#paste0(path, subfolder[j], "/", "KM_plot_train_with_PI.jpeg")
#jpeg(file="KM_plot.jpeg", units="in", width=10, height=10, res=300)
pp <-  ggsurvplot(fit_tr, data=tr_PI, pval=TRUE,
                  risk.table=TRUE, tables.height = 0.3, #add risk table & height
                  xlab="Time in Months",
                  legend.labs=c( "Less than Mean of PI", "Greater than Mean of PI"), legend.title="PI Score",  
                  palette=c("dodgerblue2", "red"), 
                  risk.table.col="strata", break.time.by = 12,
                  conf.int = F, censor = TRUE,
                  surv.median.line = "hv", # Add medians survival
                  #palette = c("red","blue"),#add desired color
                  size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                  font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
pp
jpeg("Train_PI_KM_plot.jpg", units="in", width=10, height=10, res=300)
print(pp, newpage = FALSE)
dev.off()

#ROC plots

LUSC_C_tr$PI


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
jpeg(file="ROC_train.jpeg", units="in", width=10, height=10, res=300)
plot(tr_roc1$FP, tr_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),
     xlab="FP", col="red",
     ylab="TP",main= "AUC Curve for Survival Prediction")
lines(tr_roc3$FP, tr_roc3$TP, type="l", lty=2, col="blue")
lines(tr_roc5$FP, tr_roc5$TP, type="l", lty=2, col="green")
legend(0.5,0.5, legend=c( paste( "1 Year AUC = ",round(tr_roc1$AUC,2)),  paste("3 Years AUC = ",round(tr_roc3$AUC,2)),  paste("5 Years AUC = ",round(tr_roc5$AUC,2))), col =c ("red","blue", "green"), lty=c(1,2), bty="n")

abline(0,1)
dev.off()


