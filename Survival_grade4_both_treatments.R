
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
args <- commandArgs(TRUE)

setwd("/Users/kaurh8/Documents/CCGA_datasets/CGGA_mRNA_693_samples/Primary_grade4_survival/")

#setwd("/Users/kaurh8/Documents/CCGA_datasets/mRNA_693_samples/LGG/LGG/survival")

set.seed(7)

#data <- read.table(args[1], header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
data <- read.table("Primary_grade4_clin_data_STS_MTS_LTS_both_trts_with_exp_data.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
dim(data)
head(Clin_data)
head(data[1:25])
Clin_data <- data[1:24]


data1 <- subset(data, OS_month!="NA")
dim(data1)

# create a file to store results
write.table(cbind("ID","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="Survival_results_for_genes.txt",row.names=F,col.names=F,sep = '\t');

#Here features to compute survival start from 16th column onwards, which is transcriptomics features

for(i in seq(from=25, to=length(data1), by=1))
{
  surv_object2 <- Surv(time = data1$OS_month, event = data1$Death_Status)
  
  #survival analysis: fits cox ph model to find HR for mean cut
  fit1 <- survfit(surv_object2 ~ (data1[,i])>(mean(data1[1,i])), data=data1);
  summary(fit1);
  fit1.coxph <- coxph(surv_object2 ~ (data1[,i])>(mean(data1[1,i])), data = data1)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  
  write.table(cbind(colnames(data1[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
              file="Survival_results_for_genes.txt",row.names=F,col.names=F,sep = '\t',append = T);#output file
  
}


### Significant results 

write.table(cbind("ID","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="Significant_Survival_results_for_genes.txt",row.names=F,col.names=F,sep = '\t');


for(i in seq(from=25, to=length(data), by=1))
{
  surv_object2 <- Surv(time = data1$OS_month, event = data1$Death_Status)
  
  #survival analysis: fits cox ph model to find HR for mean cut
  fit1 <- survfit(surv_object2 ~ (data1[,i])>(mean(data1[1,i])), data=data1);
  summary(fit1);
  fit1.coxph <- coxph(surv_object2 ~ (data1[,i])>(mean(data1[1,i])), data = data1)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  
  #check whether the pvalue is significant (< 0.05) or not
  if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
  {write.table(cbind(colnames(data1[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file="Significant_Survival_results_for_genes.txt",row.names=F,col.names=F,sep = '\t',append = T);#output file
  }
}




##############
# create survival object
surv_object_p <- Surv(time = data1$OS_month, event = data1$Death_Status)
head(data1[1:13])

#Cox multivariate model
cox_multivariate <- coxph(surv_object_p ~  Gender + Grade +Group +  as.numeric(Age) + IDH_mutation_status + Codel_1p19q_status + MGMTp_meth_status  ,  data=data1 )
summary(cox_multivariate )
jpeg(file="multivariate.jpeg", units="in", width=10, height=10, res=300)
#ggforest(cox_multivariate,data=data1 )
ggforest(cox_multivariate,data=data1)
dev.off()

covariates <- c( "Gender" , "Grade", "Group",  "Age", "IDH_mutation_status", "Codel_1p19q_status", "MGMTp_meth_status")
univ_formulas <- sapply(covariates,  function(x) as.formula(paste('surv_object_p ~', x)))

univ_models <- lapply(univ_formulas, function(x){coxph(x, data = data1 )})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #p.value<-signif(x$wald["pvalue"], digits=2) 
                         #wald<-signif(x$wald["test"], digits=2)
                         p.value<-signif(x$logtest["pvalue"], digits=2)
                         logrank<-signif(x$logtest["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         CI1 <- round(as.numeric(signif(x$concordance[1])),2);#concordance
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, CI1, logrank, p.value)
                         #res<-c(beta, HR,  wald, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "CI", "logrank",  "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })

univ_results
univ_results
res_uni <- as.matrix(univ_results)
res_uni

write.table(res_uni , file = "Univariate_results_traing_data.txt", sep="\t", quote=F, row.names = T)


# Gender based
fit_gender_p = survfit(surv_object_p ~ Gender, data=data1)

?ggsurvplot
# we produce a Kaplan Meier plot
jpeg(file="Gender.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_gender_p, data=data1, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.mean.line = "hv", # Add means survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()
fit.coxph_gender_p <- coxph(surv_object_p ~ Gender, data=data1)
coef(summary(fit.coxph_gender_p))





# Class - R/P based
fit_subtype_p = survfit(surv_object_p ~ Group, data=data1)

# we produce a Kaplan Meier plot
jpeg(file="Subtype.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_subtype_p, data=data1, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()
fit.coxph_subtype_p <- coxph(surv_object_p ~ Group, data=data1)
coef(summary(fit.coxph_subtype_p))



#remove samples with no IDH mutation information
GBM_data_idh_p<-subset(data1, IDH_mutation_status!="NA")

#Create survival object
surv_object_idh_p <- Surv(time = GBM_data_idh_p$OS_month, event = GBM_data_idh_p$Death_Status)

#fit model
fit_idh_p = survfit(surv_object_idh_p ~ IDH_mutation_status, data=GBM_data_idh_p)

# Kaplan Meier plot
ggsurvplot(fit_idh_p, data=GBM_data_idh_p, pval=TRUE, risk.table=TRUE, xlab="Time in months", risk.table.col="strata")
jpeg(file="IDH1mut.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_idh_p, data=GBM_data_idh_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()

fit.coxph_idh_p <- coxph(surv_object_idh_p ~ IDH_mutation_status, data=GBM_data_idh_p)
coef(summary(fit.coxph_idh_p))

#remove samples with no MGMT information
GBM_data_mgmt_p<-subset(data1, MGMTp_meth_status!="NA")

#Create survival object
surv_object_mgmt_p <- Surv(time = as.numeric(GBM_data_mgmt_p$OS_month), event = GBM_data_mgmt_p$Death_Status)

#fit model
fit_mgmt_p = survfit(surv_object_mgmt_p ~ MGMTp_meth_status, data=GBM_data_mgmt_p)

# Kaplan Meier plot
ggsurvplot(fit_mgmt_p, data=GBM_data_mgmt_p, pval=TRUE, risk.table=TRUE, xlab="Time in days", risk.table.col="strata")

jpeg(file="mgmt.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_mgmt_p, data=GBM_data_mgmt_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()


fit.coxph_mgmt_p <- coxph(surv_object_mgmt_p ~ MGMTp_meth_status, data=GBM_data_mgmt_p)
coef(summary(fit.coxph_mgmt_p))



#remove samples with no codeletion mutation information
GBM_data_codel_p<-subset(data1, Codel_1p19q_status!="NA")

#Create survival object
surv_object_codel_p <- Surv(time = as.numeric(GBM_data_codel_p$OS_month), event = GBM_data_codel_p$Death_Status)

#fit model
fit_codel_p = survfit(surv_object_codel_p ~ Codel_1p19q_status, data=GBM_data_codel_p)

# Kaplan Meier plot
ggsurvplot(fit_codel_p, data=GBM_data_codel_p, pval=TRUE, risk.table=TRUE, xlab="Time in months", risk.table.col="strata")
jpeg(file="codelmut.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_codel_p, data=GBM_data_codel_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()


fit.coxph_codel_p <- coxph(surv_object_codel_p ~ Codel_1p19q_status, data=GBM_data_codel_p)
coef(summary(fit.coxph_codel_p))


#remove samples with no Survival class information
GBM_data_survC_p<-subset(data1, Survival_class!="NA")

#Create survival object
surv_object_survC_p <- Surv(time = as.numeric(GBM_data_survC_p$OS_month), event = GBM_data_survC_p$Death_Status)

#fit model
fit_survC_p = survfit(surv_object_survC_p ~ Survival_class, data=GBM_data_survC_p)

#summary(fit_survC_p)

# Kaplan Meier plot
ggsurvplot(fit_survC_p, data=GBM_data_survC_p, pval=TRUE, risk.table=TRUE, xlab="Time in days", risk.table.col="strata")

jpeg(file="Survival_class.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_survC_p, data=GBM_data_survC_p, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.3, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
dev.off()


fit.coxph_survC_p <- coxph(surv_object_survC_p ~ Survival_class, data=GBM_data_survC_p)
coef(summary(fit.coxph_survC_p))


#GBM_p_surv_data_age <- subset(GBM_p_surv_data, Age_group!="pediatric")
GBM_p_surv_data_age = subset(data1, Age!="NA")
#Create survival object
surv_object_p_age <- Surv(time = as.numeric(GBM_p_surv_data_age$OS_month), event = GBM_p_surv_data_age$Death_Status)

#fit model
fit_age = survfit(surv_object_p_age ~ Age_group, data=GBM_p_surv_data_age)

# Kaplan Meier plot
ggsurvplot(fit_age , data=GBM_p_surv_data_age, pval=TRUE, risk.table=TRUE, tables.height = 0.4, xlab="Time in months", risk.table.col="strata")

jpeg(file="age_groups.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_age, data=GBM_p_surv_data_age, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.4, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))

dev.off()
fit.coxph_surv_age_p <- coxph(surv_object_p_age ~ Age_group, data=GBM_p_surv_data_age)
coef(summary(fit.coxph_surv_age_p))


surv_object_p_age <- Surv(time = as.numeric(GBM_p_surv_data_age$OS_month), event = GBM_p_surv_data_age$Death_Status)

#fit model
fit_age = survfit(surv_object_p_age ~ Age_group, data=GBM_p_surv_data_age)

# Kaplan Meier plot
ggsurvplot(fit_age , data=GBM_p_surv_data_age, pval=TRUE, risk.table=TRUE, tables.height = 0.4, xlab="Time in months", risk.table.col="strata")

jpeg(file="age_groups.jpeg", units="in", width=10, height=10, res=300)
ggsurvplot(fit_age, data=GBM_p_surv_data_age, pval=TRUE, 
           risk.table=TRUE, tables.height = 0.4, #add risk table & height
           xlab="Time in months", 
           risk.table.col="strata", break.time.by = 12,
           conf.int = F, censor = TRUE,
           surv.median.line = "hv", # Add medians survival
           #palette = c("red","blue"),#add desired color
           size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"), 
           font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))

dev.off()
fit.coxph_surv_age_p <- coxph(surv_object_p_age ~ Age, data=GBM_p_surv_data_age)
coef(summary(fit.coxph_surv_age_p))

clin_factors_Coeff <- rbind(coef(summary(fit.coxph_gender_p)), coef(summary(fit.coxph_subtype_p)),  coef(summary(fit.coxph_surv_age_p)), coef(summary(fit.coxph_idh_p)), coef(summary(fit.coxph_codel_p)),  coef(summary(fit.coxph_mgmt_p)))
clin_factors_Coeff 


#Write into files
write.table(clin_factors_Coeff, file = "clin_factors_Coeff_train_COXPH.txt", sep="\t", quote=F, row.names = T)




######## cox modellasso model ###


#data1 <- read.table("282_Sig_genes_Quantile_mat_data_with_clin.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
data1 <- read.table("171_Sig_genes_Quantile_mat_data_with_clin.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

surv_object2 <- Surv(time = data1$OS_month, event = data1$Death_Status)
cvfit1 <- cv.glmnet(as.matrix(data1[25:ncol(data1)]), 
                    surv_object2, # create survival object from the data
                    family = "cox", # specify Cox PH model
                    type.measure = "C", 
                    nfolds = 5, 
                    alpha = 1, # lasso: alpha = 1; ridge: alpha=0
                    maxit = 1000)

lambda_min <- cvfit1$lambda.min
lambda_min

plot(cvfit1)


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
key_variables


write.table(key_variables,file="Lasso_key_variables.txt",row.names=T,col.names=T,sep = '\t', quote = F)



PI_model <- AKR1C3 * (-0.148) + ATOH8 *(-0.127) + CP * (0.129) + CRLF1 * (-0.157) + DRAXIN * (0.325) + GLIS1 * (-0.48) + H19 *(0.0186) + HIST1H2AG * (0.2258) +` RP11-189B4.6`* (0.3927) + `RP11-300M24.1` * (0.1277) + TNFSF14 * (0.1389)

data2 <-read.table("11_lasso_Sig_genes_Quantile_mat_data_with_clin_copy.txt",header =TRUE, sep = "\t", row.names=1,  check.names = FALSE)
######## Load file containing results (significant genes with p-value <0.05) from univariate analysis with beta coefficient  ##############################
dim(data2)
sel_features_results<-read.table("11_key_features.txt",header =TRUE, sep = "\t", check.names = FALSE)
head(sel_features_results)
dim(sel_features_results)
sel_ftr_surv <- as.data.frame(sel_features_results[,1])
names(sel_ftr_surv ) <- c("ID")
head(sel_ftr_surv )
dim(sel_ftr_surv )

#### prepare training, test and external validation data with selected features having significant value in univariate analysis ##########
sel_train <- as.data.frame(data2[,colnames(data2) %in% c(sel_ftr_surv$ID), ])
head(sel_train,2)
dim(sel_train)

data2$OS_month
######### Make final files with the selected features (genes here)  and combine survival information ###########################################################
train_feature_mat<- cbind(data2["OS_month"],data2["Death_Status"],sel_train)
head(train_feature_mat,2)
dim(train_feature_mat)

############ remove where OS.time=NA ############
train_feature_mat1<-subset(train_feature_mat,OS_month!="NA")
train_feature_mat1<-subset(train_feature_mat1,OS_month!=0)


train_feature_mat2<- cbind(data2["OS_month"],data2["Group"],sel_train)
# save files with selected genes & survival information #########
write.table(cbind("ID"=rownames(train_feature_mat1), train_feature_mat1),file="sel_train.txt",sep="\t",quote=F, row.names=F)

#Create prognostic index for pathway
LUSC_C_tr=train_feature_mat1
tr <- LUSC_C_tr[3:ncol(LUSC_C_tr)]

E= length(tr)
E

PI_tr=0
for(i in seq(from=1, to=E ,by=1))
{
  PI_tr= PI_tr+((tr[,i])*(sel_features_results[i,2]))
}


LUSC_C_tr$PI<-PI_tr
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
                  risk.table.col="strata", break.time.by = 12,
                  conf.int = F, censor = TRUE,
                  surv.median.line = "hv", # Add medians survival
                  #palette = c("red","blue"),#add desired color
                  size=1.5, font.tickslab = c(12,  "black"), font.y = c(14,  "black"),
                  font.x = c(14, "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"))
pp
jpeg("KM_plot.jpg", units="in", width=10, height=10, res=300)
print(pp, newpage = FALSE)
dev.off()



# multivariate
tr_PI_data_with_clin.txt