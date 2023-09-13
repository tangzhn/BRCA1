######Video source: http://ke.biowolf.cn
######Homepage: http://www.biowolf.cn/
######Wechat Public：biowolf_cn
######E-mail：biowolf@foxmail.com
######Wechat: 18520221056

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSEABase")

#setwd("C:\\Users\\biowolf\\Desktop\\ssGSEA\\07.ssGSEA")          #设置工作目录
inputFile="uniq.symbol.txt"                                         #输入文件
gmtFile="immune.gmt"                                           #GMT文件

#引用包
library(GSVA)
library(limma)
library(GSEABase)

#读取输入文件，并对输入文件处理
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt[1:3,1:3]
#rt=as.matrix(rt)
#rownames(rt)=rt[,1]
#exp=rt[,2:ncol(rt)]
#dimnames=list(rownames(exp),colnames(exp))
#mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
#mat=avereps(mat)
#mat=mat[rowMeans(mat)>0,]


rownames(rt)=rt[,1]
rt <- rt[,-1]
mat <- as.matrix(rt)

geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())

#ssgsea分析
ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#定义ssGSEA score矫正函数
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
#对ssGSEA score进行矫正
ssgseaOut=normalize(ssgseaScore)
ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
write.table(ssgseaOut,file="ssgseaOut.txt",sep="\t",quote=F,col.names=F)

######Video source: http://ke.biowolf.cn
######Homepage: http://www.biowolf.cn/
######Wechat Public：biowolf_cn
######E-mail：biowolf@foxmail.com
######Wechat: 18520221056#    箱线图  ########
##hlaexp包含：样本信息列，分组信息列，基因表达量信息列（多列）
#hlaexp <- texp_sort[,c('Cluster',hla)]
install.packages("tidyr")
install.packages("ggplot2")
install.packages("tidyverse")
library(tidyverse)
library(tidyr)
library(ggplot2)
require(ggpubr)
require(ggsci)
require(cowplot)

#######################准备画图需要的文件格式 ##################
hlaexp<-read.table("ssgseaOut.txt",sep="\t",header=T,check.names = F)
hlaexp <- as.data.frame(t(hlaexp))
colnames(hlaexp) <- hlaexp[1,]
hlaexp <- hlaexp[-1,]

hlaexp$group=c(rep("Normal",113),rep("Tumor",1091)) 
hlaexp <- hlaexp[,c(ncol(hlaexp),1:(ncol(hlaexp)-1))]
hlaexp <-cbind(id=rownames(hlaexp),hlaexp)
write.table(hlaexp,"hlaexp.txt",sep="\t",quote = F,row.names = F)

#################### 开始画图########################
hlaexp<-read.table("hlaexp.txt",sep="\t",header=T,check.names = F)
#hlaexp <- hlaexp[,-3]
rownames(hlaexp) <- hlaexp[,1]
#rownames(hlaexp) <- hlaexp$ID 
hlaexp=hlaexp[,-1]

boxdata <- gather(hlaexp,HLA,exp,2:ncol(hlaexp))
#boxdata$HLA=gsub("\\.","\\-",boxdata$HLA)
#gather之后的boxdata包含三列信息，表达量信息列，分组信息列，基因名称信息列
boxdata[,'exp'] <- sapply(boxdata[,'exp'],function(x){log2(x+1)})

#group分组列变成factor,boxplot按照L,H排列
boxdata$group <- factor(boxdata$group,levels = c('Normal','Tumor'))

#boxdata <- as.data.frame(apply(boxdata$group,2,function(x) paste('x','risk',sep='-')))


# png('HLA.boxplot.png',width = 650,height = 350)

pdf('Immune.pdf',width =9,height =6)
ggboxplot(boxdata,x='HLA',y='exp',color ='group',
          ylab = 'Relative Expression',
          xlab ='',palette = "d3",  #nejm #jco #jama#npg#d3#aaas
          add ='jitter',add.params = list(size=0.1,shape=16),   #add ='jitter'
          width=0.5)+
  rotate_x_text(45)+
  #coord_flip()+   #颠倒X轴和Y轴
  theme(axis.text=element_text(size=8),
        axis.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  stat_compare_means(size = 3,aes(group=group),
                     label = 'p.signif',method = 'wilcox.test',hide.ns = F)  ##检验方法可以改
dev.off()
