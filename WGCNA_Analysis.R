### WGCNA Analysis

library(WGCNA)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(mat, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
net = blockwiseModules(mat, power = 6,
TOMType = "unsigned", minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = "TOM",
verbose = 3)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
moduleinfo=cbind(moduleColors,row.names(mat))
MEs = net$MEs

###Gene module eigenvalue-three replicated DE-miR associations

MEsmiR <- read.delim("ME_miRNA_assocInput.txt")
# outcome
out_start=1
out_end= 20
out_nvar=out_end-out_start+1
out_variable=rep(NA, out_nvar)
out_beta=rep(NA, out_nvar)
out_OR = rep(NA, out_nvar)
out_se = rep(NA, out_nvar)
out_pvalue=rep(NA, out_nvar)
# exposure
exp_start=22
exp_end=24
exp_nvar=exp_end-exp_start+1
exp_variable=rep(NA, exp_nvar)
exp_beta=rep(NA, exp_nvar)
exp_OR = rep(NA, out_nvar)
exp_se = rep(NA, out_nvar)
exp_pvalue=rep(NA, exp_nvar)
number=1

for (i in out_start:out_end){
    outcome = colnames(MEsmiR)[i]
    for (j in exp_start:exp_end){
        exposure = colnames(MEsmiR)[j]
        model <- lmer(get(outcome) ~ get(exposure) + (1 | PatientID) + SEX + AGE + RACE
                      ,
                      na.action = na.exclude,
                      data=MEsmiR)
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
LMM_ME_miRNA=cbind(expo,expo_p.adjust,out,out_p.adjust)
library(reshape2)
OR=acast(LMM_ME_miRNA,out_variable~exp_variable,value.var="exp_OR")
adjP=acast(LMM_ME_miRNA,out_variable~exp_variable,value.var="expo_p.adjust")


textMatrix = paste(signif(OR, 2), "\n(",signif(adjP, 1), ")", sep = "");
dim(textMatrix) = dim(OR)
par(mar = c(6, 20, 3, 1));
colfunc <- colorRampPalette(c("indianred", "white"))
labeledHeatmap(Matrix = adjP,
               xLabels = colnames(adjP),
               yLabels = row.names(adjP),
               ySymbols = row.names(adjP),
               colorLabels = FALSE,
               colors = colfunc(10),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.9,
               zlim = c(0,1),
               main = paste("Module-miRNA relationships"),cex.lab.x = 0.8,
               cex.lab.y = 0.8,font.lab.x=2,font.lab.y=2)
