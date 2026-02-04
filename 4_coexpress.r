##计算籽粒基因表达相关性
data<- read.table("rna_new_residuals.txt",header = T,check.names = F,row.names = 1)
data <- as.matrix(data)
data = data[apply(data, 1, function(x) sd(x)!=0),] 
data = data[,apply(data, 2, function(x) sd(x)!=0)] 
library(Hmisc)
options(digits = 2)	
cor_dat2=rcorr(data,type="pearson")
alpha <- 0.05
cor_dat2$r[cor_dat2$P > alpha] <- NA
cor_dat2$r=round(cor_dat2$r,2)  # 相关系数
cor_dat2$P=round(cor_dat2$P,2)  # P值
write.table(cor_dat2$r,file = "kernel-6day_corexpress_r.txt",row.names = T,sep = "\t",quote = T)
write.table(cor_dat2$P,file = "kernel_6day_corexpress_p.txt",row.names = T,sep = "\t",quote = T)
