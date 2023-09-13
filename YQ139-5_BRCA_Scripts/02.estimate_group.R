

library(ggplot2)
library(ggpubr)
library(magrittr)
library(ggsignif)

#inputFile="clinical_treated.txt"                                      
#rt=read.table(inputFile,sep="\t",header=T,check.names=F)

#inputFile="临床和风险因子综合表.csv"                                      
#rt=read.table(inputFile,sep=",",header=T,check.names=F)
rt=read.table("scoersClinical.txt",sep="\t",header=T,check.names=F,row.names = 1)
colnames(rt)

###igv,jama
#4组
rt <- rt[order(rt$T),]
pdf("T_riskScore.box.pdf",height=6,width=6)
ggboxplot(rt, x= 'T', y='ImmuneScore',
          ylab = "ImmuneScore",xlab = "", color = 'T', palette = "aaas", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= T))
dev.off()

################    N重新画
rt=read.table("new_N.txt",sep="\t",header=T,check.names=F,row.names = 1)
unique(rt$N)
colnames(rt)

rt <- rt[order(rt$N),]
pdf("N_new_riskScore.box.pdf",height=6,width=6)
ggboxplot(rt, x= 'N', y='ImmuneScore',
          ylab = "ImmuneScore",xlab = "", color = 'N', palette = "aaas", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= N))
dev.off()



t <- rt[order(rt$neoplasm_histologic_grade),]   #为T升序排列
pdf("grade_riskScore.box.pdf",height=6,width=6)
ggboxplot(t, x= 'neoplasm_histologic_grade', y='riskScore',
          ylab = "riskScore", xlab = "",color = 'neoplasm_histologic_grade',  palette = "futurama", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= neoplasm_histologic_grade))
dev.off()


inputFile="临床和风险因子综合表.csv"                                      
rt=read.table(inputFile,sep=",",header=T,check.names=F)
#rt <- rt[order(rt$tumor_residual_disease),];unique(rt$tumor_residual_disease)
#rt1 <- subset(rt,N!="N3")   ###删除N列是N3的行   去除数据中有NA的行：ge <- na.omit(ge)
#unique(rt1$N)
pdf("tumor_residual_disease.box.pdf",height=6,width=6)
ggboxplot(rt, x= 'tumor_residual_disease', y='riskScore',
          ylab = "riskScore", xlab="",color = 'tumor_residual_disease', palette = "aaas", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= tumor_residual_disease))
dev.off()

#3组



##2组
rt <- rt[order(rt$race),];unique(rt$race)
rt1 <- subset(rt,race!="american indian or alaska native")

pdf("race.box.pdf",height=6,width=6)
ggboxplot(rt1, x= 'race', y='riskScore',
          #ylab = "FYN expression", color = 'M', palette = "igv", merge = "flip", add="jitter")+
          ylab = "riskScore", xlab="",color = 'race', palette = "aaas", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= race))
dev.off()



rt <- rt[order(rt$gender),]
pdf("gender.box.pdf",height=6,width=6)
ggboxplot(rt, x= 'gender', y='riskScore',
          ylab = "riskScore", color = 'gender', palette = "startrek", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= gender))
dev.off()

#年龄
inputFile="临床和风险因子综合表.csv"                                      
rt=read.table(inputFile,sep=",",header=T,check.names=F)
rt$age <- ifelse(rt$age>65,'>65','<=65')
rt1 <- rt[order(rt$age),]   #为age升序排列
pdf("age.box.pdf",height=6,width=6)
ggboxplot(rt1, x= 'age', y='riskScore',
          ylab = "risk Score", xlab="",color = 'age', palette = "aaas", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= age))
dev.off()
