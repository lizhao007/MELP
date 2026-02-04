library(rrBLUP)
#random population of 200 lines with 1000 markers
M <- matrix(rep(0,200*1000),1000,200)
for (i in 1:200) {
  M[,i] <- ifelse(runif(1000)<0.5,-1,1)
}
colnames(M) <- 1:200
geno <- data.frame(marker=1:1000,chrom=rep(1,1000),pos=1:1000,M,check.names=FALSE)
QTL <- 100*(1:5) #pick 5 QTL
u <- rep(0,1000) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(M,u))
h2 <- 0.5
y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))
pheno <- data.frame(line=1:200,y=y)
scores <- GWAS(pheno,geno,plot=FALSE)


#????
setwd("D:/10 Network-driven GWAS/twas")
rm(list = ls())
library(rrBLUP)
library(qqman)
#load the Markers and Phenotypes
Markers <- as.matrix(read.table(file="ExpMatrixWW_V4-test1.txt"), header=F)
dim(Markers)
Markers[1:10,1:10]
kin <- A.mat(
  Markers, #????Ϊ{-1,0,1}?ı??Ǿ??󣬿ɰ???С???㣨????????????????NA
  min.MAF=NULL, #??С??λ????Ƶ??
  max.missing=NULL, #????ȱʧֵ????
  impute.method="mean", #???䷽????mean/EM
  tol=0.02, #EM?㷨????��??׼
  n.core=1, #EM?㷨?Ĳ??к?????????unix??
  shrink=FALSE, #?????��?
  return.imputed=FALSE #???????????Ǿ???
)
#?????궼??NA??ֻʶ??-1,0,1
trait <- read.table("eucidean_mph-125.txt",header = T)
TPM <- read.table("ExpMatrixWW_V4-125.txt",header = T,row.names = 1)
kin <- read.table("125-kin.txt",header = T,row.names = 1)
#rownames(kin) <- trait[,1]
#head(trait)
#colnames(kin) <- trait[,1]
kin <- as.matrix(kin)
pos <- read.table("exgeneid.bed")
pos <- subset(pos,pos$V4 %in% rownames(TPM))
colnames(TPM) <- rownames(kin)
exp <- data.frame(maker=rownames(TPM),chrom=pos[,1],pos=pos[,2],TPM)
######????????
head(trait[,c(1,2)])
data <- GWAS(trait[,c(1,2)], exp, plot = F, P3D = F, K=kin)
data <- GWAS(trait[,c(1,2)], exp, plot = T, P3D = T) #??֧???ṩ??K?ļ?????֧??P3D=F
data1 <- GWAS(trait[,c(1,2)], exp, plot = T, P3D = T, K=kin)
head(data)
write.table(data,file = "TWAS_M04.txt",row.names = F,sep = "\t",quote = F)
qq(10^(-data$MPH))


#twas????
setwd("D:/10 Network-driven GWAS/twas")
rm(list = ls())
library(rrBLUP)
#trait <- read.table("eucidean_mph-125.txt",header = T)
trait <- read.table("gwas_ph-125.txt",header = T)
TPM <- read.table("ExpMatrixWW_V4-125.txt",header = T,row.names = 1)
kin <- read.table("kinship-125.txt",header = F)
rownames(kin) <- trait[,1]
colnames(kin) <- trait[,1]
kin <- as.matrix(kin)
pos <- read.table("exgeneid.bed")
pos <- subset(pos,pos$V4 %in% rownames(TPM))
exp <- data.frame(maker=rownames(TPM),chrom=pos[,1],pos=pos[,2],TPM)
# ͳ??exp??????ÿ??Ϊ0??Ԫ?ظ???
n0 <- apply(exp == 0, 1, sum)
# ͳ??0Ԫ?ش???100??80%???????к?
i0 <- which(n0 > 100) 
# ɾ??0Ԫ?ش???100??????
exp <- exp[-i0, ]
data <- GWAS(trait, exp, plot = T, P3D = T, K=kin) #??֧???ṩ??K?ļ?????֧??P3D=F
#write.table(data,file = "TWAS_125-filter.txt",row.names = F,sep = "\t",quote = F)
write.table(data,file = "TWAS_125-filter-gwas-ph.txt",row.names = F,sep = "\t",quote = F)




#twas????--200line-100grainweight
setwd("D:/10 Network-driven GWAS/kernel/twas")
rm(list = ls())
library(rrBLUP)
trait <- read.table("100grainweight-190line.txt",header = T, check.names = F)
TPM <- read.table("rna_new_residuals-190.txt",header = T, check.names = F,row.names = 1)
kin <- read.table("190-kinship-for-rrblup.txt",header = F)
rownames(kin) <- trait[,1]
colnames(kin) <- trait[,1]
kin <- as.matrix(kin)
pos <- read.table("rna-gene.bed")
pos <- subset(pos,pos$V1 %in% rownames(TPM))
exp <- data.frame(maker=rownames(TPM), check.names = F,chrom=pos[,2],pos=pos[,3],TPM)
data <- GWAS(trait, exp, plot = T, P3D = T, K=kin) #??֧???ṩ??K?ļ?????֧??P3D=F
write.table(data,file = "190-100grainweight-for-rrblup.txt",row.names = F,sep = "\t",quote = F)




#PD-雄穗分支数-twas分析
setwd("D:/10 Network-driven GWAS/twas")
rm(list = ls())
library(rrBLUP)
#trait <- read.table("eucidean_mph-125.txt",header = T)
trait <- read.table("gwas-tasselbranchnumber-125.txt",header = T)
TPM <- read.table("ExpMatrixWW_V4-125.txt",header = T,row.names = 1)
kin <- read.table("kinship-125.txt",header = F)
rownames(kin) <- trait[,1]
colnames(kin) <- trait[,1]
kin <- as.matrix(kin)
pos <- read.table("exgeneid.bed")
pos <- subset(pos,pos$V4 %in% rownames(TPM))
exp <- data.frame(maker=rownames(TPM),chrom=pos[,1],pos=pos[,2],TPM)
# 统计exp矩阵中每行为0的元素个数
n0 <- apply(exp == 0, 1, sum)
# 统计0元素大于100（80%）个的行号
i0 <- which(n0 > 100) 
# 删除0元素大于100个的行
exp <- exp[-i0, ]
data <- GWAS(trait, exp, plot = T, P3D = T, K=kin) #不支持提供的K文件，不支持P3D=F
#write.table(data,file = "TWAS_125-filter.txt",row.names = F,sep = "\t",quote = F)
write.table(data,file = "TWAS_125-filter-gwas-tasselbranchnumber.txt",row.names = F,sep = "\t",quote = F)



