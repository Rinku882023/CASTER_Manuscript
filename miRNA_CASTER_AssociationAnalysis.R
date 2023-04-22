### Investigating BATCH EFFECT 

library(gPCA)
rawCount <- read.delim("rawCount.txt", row.names=1)
target <- read.delim("target.txt")
data=as.matrix(t(rawCount))
batch= as.factor(target$BATCH_3_GROUPING)
out<-gPCA.batchdetect(x=data,batch=batch,center=FALSE,nperm=1000)
out$delta ; out$p.val
gDist(out)

### DeSeq2 Normalization and shifted logarithm transformation (log2(n + 1))

library(DESeq2)
rc=as.matrix(rawCount)
all(target$S_SUBJECTID %in% colnames(rc))
all(target$S_SUBJECTID %in% colnames(rc))
all(target$S_SUBJECTID ==colnames(rc))
target$CASTER=as.factor(target$CASTER_a)
target$gender=as.factor(target$gender) 
dds <- DESeqDataSetFromMatrix(countData=round(rc), colData = target,design = ~ CASTER+age+gender)
dds1 <- estimateSizeFactors(dds)
norm <- data.frame(log2(counts(dds1, normalized=TRUE)[, 1:580] + 1))
