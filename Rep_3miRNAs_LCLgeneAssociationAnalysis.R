#### Association Analysis between three replicated miRNAs and genes

DEG_3miR_expData <- read.delim("combined.txt", row.names=1)
library(lme4)
# outcome (Gene expression columns)
out_start=3
out_end= 4820
out_nvar=out_end-out_start+1
out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_OR = rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)
# exposure (3 miRNA expression columns) 
exp_start=4821
exp_end=4823
exp_nvar=exp_end-exp_start+1
exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_OR = rep(NA, out_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)
number=1
 
for (i in out_start:out_end){
    outcome = colnames(DEG_3miR_expData)[i]
    for (j in exp_start:exp_end){
        exposure = colnames(DEG_3miR_expData)[j]
        model <- lmer(get(outcome) ~ get(exposure) + (1 | PatientID) + SEX + AGE
                      
                      ,
                      na.action = na.exclude,
                      data=DEG_3miR_expData)
        Vcov <- vcov(model, useScale = FALSE)
        beta <- fixef(model)
        OR <-exp(beta)
        se <- sqrt(diag(Vcov))
        zval <- beta / se
        pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
        
        out_beta[number] = as.numeric(beta[2])
        out_OR[number] = as.numeric(OR[2])
        out_se[number] = as.numeric(se[2])
        out_pvalue[number] = as.numeric(pval[2])
        out_variable[number] = outcome
        number = number + 1
        
        exp_beta[number] = as.numeric(beta[2])
        exp_OR[number] = as.numeric(OR[2])
        exp_se[number] = as.numeric(se[2])
        exp_pvalue[number] = as.numeric(pval[2])
        exp_variable[number] = exposure
        number = number + 1
    }
}
outcome = data.frame(out_variable, out_beta,out_OR, out_se, out_pvalue)
exposure = data.frame(exp_variable, exp_beta,exp_OR, exp_se, exp_pvalue)
expo = na.omit(exposure)
out = na.omit(outcome)
expo_p.adjust=p.adjust(expo$exp_pvalue,method="BH")
out_p.adjust=p.adjust(out$out_pvalue,method="BH")
LMM_mRNA_miRNA=cbind(expo,expo_p.adjust,out,out_p.adjust) 
write.csv(LMM_mRNA_miRNA,file="Result_miRNA_geneAssociation.csv")
