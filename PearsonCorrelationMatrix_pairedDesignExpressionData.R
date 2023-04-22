### Pearson Correlation matrix generation for paired designed expression data

DEX<- read.delim("DEX.txt")
SHAM<- read.delim("SHAM.txt")
start=1
end=  4818
PairedCorr  <- function(Xi,Xj,Yi,Yj) 
{
    p <- length(Xi)/(length(Xi) + length(Yi))
    covX <- cov(Xi,Xj, method="pearson")
    covY <- cov(Yi, Yj, method="pearson")
    d1 <- mean(Yi)- mean(Xi)
    d2 <- mean(Yj) -mean(Xj)
    varxi <- var(Xi)
    varyi <- var(Yi)
    varxj <- var(Xj) 
    varyj <- var(Yj)
    corrZiZj = ((1-p)*covX + p*covY + p*(1-p)*d1*d2)/ sqrt(((1-p)*varxi + p*varyi + p*(1-p)*d1^2)*((1-p)*varxj + p*varyj + p*(1-p)*d2^2))
    return(corrZiZj)
}
mat <-matrix(,nrow =4818,ncol =4818)

for (i in start: end)
{
    
    for (j in start: end)
    {
        dex1 = DEX[,i]
        dex2 = DEX[,j]
        sham1 = SHAM[,i]
        sham2 = SHAM[,j]
        mat[i,j] = PairedCorr(sham1,sham2,dex1,dex2) 
        
    }
}
row.names(mat)=row.names(DEX)
colnames(mat)=row.names(DEX)
