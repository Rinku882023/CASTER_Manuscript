### Table 1. Baseline Epidemiologic and Clinical Characteristics of the cohort data
library(table1)
library(htmlTable)
library(knitr)
library(rmarkdown)
library(MatchIt)
library(readxl)

data$CASTER_binary <- factor(data$CASTER_binary,levels=c(0,1),labels=c("Poor Responder",
                                                                   "Good Responder"))
data$CASTER_binary<-factor(data$CASTER_binary,levels=c(levels(data$CASTER_binary),"P-value"))


data$gender <- factor(data$gender, levels=c(1,2),labels=c("Male","Female"))

label(data$gender)<- "Gender"
label(data$age)<- "Age"
label(data$htcm)<- "Height"
label(data$Weight_in_kg)<-"Weight"
label(data$BMI)<- "BMI"
label(data$G22D_Bursts)<- "Bursts"
label(data$G20D_EDHOSP)<- "ED (emergency department) visits"
label(data$pctpred_fev1_pre_BD_IMPUTE)<- "% predicted Pre-BD FEV1 [(PreFEV/FEVpred) *100]"
label(data$LNPC20_IMPUTE)<- "Airway hyper-responsiveness"
label(data$bdrpred_IMPUTE)<- "Bronchodilator Response as % of baseline FEV1"
label(data$log10Ige)<- "Log10 IgE"
label(data$X_25_OH_D)<- "25 Hydroxyvitamin D"
label(data$log10eos)<- "Log10  Eosinophil"

units(data$age)<- "years"
units(data$htcm)<- "cm"
units(data$Weight_in_kg)<- "kg"
units(data$X_25_OH_D)<- "ng/ml"

rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- data[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.test(y ~ data$CASTER_binary)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(data$CASTER_binary)))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}
rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

table1(~ gender + age + htcm+Weight_in_kg+BMI+G22D_Bursts+G20D_EDHOSP+pctpred_fev1_pre_BD_IMPUTE+LNPC20_IMPUTE+bdrpred_IMPUTE+log10Ige+X_25_OH_D+log10eos| CASTER_binary, droplevels=F,render=rndr, render.strat=rndr.strat,data=data,overall=F)
