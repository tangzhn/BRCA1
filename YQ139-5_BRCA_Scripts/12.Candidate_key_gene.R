
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
library(edgeR)
library(limma)
####################################去重复基因   #######
inputFile="uniq.symbol.txt"                                                  #输入文件名字
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
boxplot(rt[,2:ncol(rt)],las=2)   #查看数据分布
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]   #把样本表达量数据提取出来
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data1 <- log2(data+1)   ###标准化数据
data1 <- data
data2 <- data1[rowMeans(data1)>0,]
data <- data2
boxplot(data[,2:ncol(data)],las=2)

save(data,file="data_fpkm.Rdata")
load("data_fpkm.Rdata")
########## 2.差异分析部分   #############
#data_FPKM <- as.data.frame(data)
data_FPKM <- read.table("uniq.symbol.txt",sep="\t",header=T,check.names = F)
data_FPKM[1:3,1:3]
#data_FPKM <- read.csv("matrix_VaD_Control.csv",header=T,check.names = F)
#data_FPKM <- data_FPKM[rowMeans(data_FPKM)>0,]
rownames(data_FPKM) <- data_FPKM[,1]
data_FPKM <- data_FPKM[,-1]
#data_FPKM <- data_FPKM[,-35:-59]  #删除了25个治疗的样本
boxplot(data_FPKM[,1:ncol(data_FPKM)],las=2)    ##查看数据分布，看是否进行标准化
#data_FPKM<-log2(data_FPKM+1)   #标准化
##data_FPKM<-rpkm(data_FPKM,gene.length =2000,normalized.lib.sizes = FALSE,lib.size = NULL,log=TRUE,prior.count =2)

##### 样本不按顺序的时候

group <- colnames(data_FPKM)
annotation_col=data.frame(group=rep(c("Normal","Tumor"),c(113,1091)))
annotation_col_1 <- data.frame(t(annotation_col))
clinical <- rbind(group,annotation_col_1)

#clinical<-read.table("group.txt",sep="\t",header=F)
clinical<-t(clinical)
clinical<-data.frame(clinical)
colnames(clinical)<-c("var1","var2")
rownames(clinical)<-clinical$var1

clinical<-clinical[as.vector(colnames(data_FPKM)),]  #如果跑不通这步就#掉
group<-factor(clinical$var2)
design <- model.matrix(~0+group, data=clinical)
colnames(design) <- c(levels(group))
identical(colnames(data_FPKM),rownames(clinical))

fit<-lmFit(data_FPKM,design)
contrast.matrix <- makeContrasts("Tumor-Normal",levels=design)  ##看clinical文件的命名
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
result1<-topTable(fit2,coef=1,adjust="BH",number=dim(data_FPKM)[1])

#write.table(result1,file="DEGs_information.txt",sep="\t",quote=F,row.names=T)
write.csv(result1,file="DEGs_information.csv",quote=F,row.names=T)
save(result1,file="差异结果.Rdata")
load("差异结果.Rdata")
#a <- result1[which(result1$P.Value<0.05 & abs(result1$logFC)>1),]

#########################   logFc=2,adj.P.Val=0.05   #######################
colnames(result1)
logFc=1
p=0.05
diffSig = result1[(result1$P.Value < p & (result1$logFC>logFc | result1$logFC<(-logFc))),]
diffSig <- cbind(gene=rownames(diffSig),diffSig)
write.table(diffSig, file="diffSig_logFC1_P0.05.xls",row.names=F,sep="\t",quote=F)
#UP
diffUp = result1[(result1$P.Value < p & (result1$logFC>logFc)),]
diffUp <- cbind(gene=rownames(diffUp),diffUp)
write.table(diffUp, file="up_logFC1_p0.05.xls",row.names=F,sep="\t",quote=F)
#Down
diffDown = result1[(result1$P.Value < p & (result1$logFC<(-logFc))),]
diffDown <- cbind(gene=rownames(diffDown),diffDown)
write.table(diffDown, file="down_logFC1_p0.05.xls",row.names=F,sep="\t",quote=F)

save(data_FPKM,diffSig,clinical,file="logistic_multibox_input.Rdata")
load("logistic_multibox_input.Rdata")

###导出差异表达基因矩阵
diffExp=data_FPKM[rownames(diffSig),]
d <- cbind(id=rownames(diffExp),diffExp)
write.table(d,file="diffmRNAExp.txt",sep="\t",quote=F,row.names=F)

#####转置差异表达基因矩阵，为了制备单因素输入文件
d_t <- t(diffExp)
dd <- cbind(id=rownames(d_t),d_t)
write.table(dd,file="Degs_Exp.txt",sep="\t",quote=F,row.names=F)


######制备multiBox_input和lasso-SVM输入数据
d_t <- as.data.frame(t(diffExp))
dd <- cbind(id=rownames(d_t),d_t)

colnames(clinical) <- c("id","group")
merg <- merge(clinical,dd,by="id",all=F)
order <- merg[order(merg$group),]
write.table(order,file="lasso_SVM_input.txt",sep="\t",quote=F,row.names=F)




if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
library(edgeR)
library(limma)
########################  1.获得部分样本的表达量   ########################
whole  <- read.table("C:\\Users\\Administrator\\Desktop\\YQ139-5\\08_Tumor-vs-Normal及亚型间差异分析\\01-Tumor-vs-Normal\\uniq.symbol.txt",sep="\t",header=T,row.names = 1,check.names = F)
whole[1:3,1:3]
whole_t <- data.frame(t(whole),check.names = F)

#######  挑出肿瘤样本
whole_t$type <- ifelse(substr(rownames(whole_t),14,15)  =="01","Tumor","Normal")
table(whole_t$type)
DEimmune1 <- subset(whole_t,whole_t$type=="Tumor")

DEimmune2  <- DEimmune1[,-ncol(DEimmune1)]  
DEimmune2$id <- substr(rownames(DEimmune2),1,12)      ###命名和time一致
DEimmune2 <- DEimmune2[,c(ncol(DEimmune2),1:(ncol(DEimmune2)-1))]  ##也可tarexp <- tarexp %>%  select(id,everything())
rownames(DEimmune2)


#############  2.分组文件   ##################
sample <- read.table("three_group.txt",sep="\t",header=T,check.names = F)
sample[1:3,1:3]
colnames(sample)

sample1 <- sample[,c(1,ncol(sample))];sample1[1:3,1:3]
sample2 <- sample1[sample1$Subtype!="C3",];table(sample2$Subtype)
sample2 <- sample2[order(sample2$Subtype),]

DEimmune2[1:3,1:3]
matrix <- DEimmune2[sample2$id,]

identical(matrix$id,sample2$id)

####################3.分别提取matrix矩阵文件和group分组文件  #####################

rownames(matrix) <- matrix$id
matrix[1:3,1:3]
matrix1 <- matrix[,-1]
data_FPKM <- data.frame(t(matrix1),check.names = F)


########## 2.差异分析部分   #############
#data_FPKM <- as.data.frame(data)
#data_FPKM <- read.table("uniq.symbol.txt",sep="\t",header=T,check.names = F)
data_FPKM[1:3,1:3]
#data_FPKM <- read.csv("matrix_VaD_Control.csv",header=T,check.names = F)
#data_FPKM <- data_FPKM[rowMeans(data_FPKM)>0,]
#rownames(data_FPKM) <- data_FPKM[,1]
#data_FPKM <- data_FPKM[,-1]
#data_FPKM <- data_FPKM[,-35:-59]  #删除了25个治疗的样本
boxplot(data_FPKM[,1:ncol(data_FPKM)],las=2)    ##查看数据分布，看是否进行标准化
#data_FPKM<-log2(data_FPKM+1)   #标准化
##data_FPKM<-rpkm(data_FPKM,gene.length =2000,normalized.lib.sizes = FALSE,lib.size = NULL,log=TRUE,prior.count =2)

identical(colnames(data_FPKM),sample2$id)
##### 样本不按顺序的时候

group <- colnames(data_FPKM)
annotation_col=data.frame(group=rep(c("C1","C2"),c(424,500)))
annotation_col_1 <- data.frame(t(annotation_col))
clinical <- rbind(group,annotation_col_1)

#clinical<-read.table("group.txt",sep="\t",header=F)
clinical<-t(clinical)
clinical<-data.frame(clinical)
colnames(clinical)<-c("var1","var2")
rownames(clinical)<-clinical$var1

clinical<-clinical[as.vector(colnames(data_FPKM)),]  #如果跑不通这步就#掉
group<-factor(clinical$var2)
design <- model.matrix(~0+group, data=clinical)
colnames(design) <- c(levels(group))
identical(colnames(data_FPKM),rownames(clinical))

fit<-lmFit(data_FPKM,design)
contrast.matrix <- makeContrasts("C1-C2",levels=design)  ##看clinical文件的命名
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
result1<-topTable(fit2,coef=1,adjust="BH",number=dim(data_FPKM)[1])

#write.table(result1,file="DEGs_information.txt",sep="\t",quote=F,row.names=T)
write.csv(result1,file="DEGs_information_C1-vs-C2.csv",quote=F,row.names=T)
save(result1,file="差异结果.Rdata")
load("差异结果.Rdata")
#a <- result1[which(result1$P.Value<0.05 & abs(result1$logFC)>1),]

#########################   logFc=2,adj.P.Val=0.05   #######################
colnames(result1)
logFc=1
p=0.05
diffSig = result1[(result1$P.Value < p & (result1$logFC>logFc | result1$logFC<(-logFc))),]
diffSig <- cbind(gene=rownames(diffSig),diffSig)
write.table(diffSig, file="diffSig_logFC1_P0.05.xls",row.names=F,sep="\t",quote=F)
#UP
diffUp = result1[(result1$P.Value < p & (result1$logFC>logFc)),]
diffUp <- cbind(gene=rownames(diffUp),diffUp)
write.table(diffUp, file="up_logFC1_p0.05.xls",row.names=F,sep="\t",quote=F)
#Down
diffDown = result1[(result1$P.Value < p & (result1$logFC<(-logFc))),]
diffDown <- cbind(gene=rownames(diffDown),diffDown)
write.table(diffDown, file="down_logFC1_p0.05.xls",row.names=F,sep="\t",quote=F)

save(data_FPKM,diffSig,clinical,file="logistic_multibox_input.Rdata")
load("logistic_multibox_input.Rdata")

###导出差异表达基因矩阵
diffExp=data_FPKM[rownames(diffSig),]
d <- cbind(id=rownames(diffExp),diffExp)
write.table(d,file="diffmRNAExp.txt",sep="\t",quote=F,row.names=F)

#####转置差异表达基因矩阵，为了制备单因素输入文件
d_t <- t(diffExp)
dd <- cbind(id=rownames(d_t),d_t)
write.table(dd,file="Degs_Exp.txt",sep="\t",quote=F,row.names=F)


######制备multiBox_input和lasso-SVM输入数据
d_t <- as.data.frame(t(diffExp))
dd <- cbind(id=rownames(d_t),d_t)

colnames(clinical) <- c("id","group")
merg <- merge(clinical,dd,by="id",all=F)
order <- merg[order(merg$group),]
write.table(order,file="lasso_SVM_input.txt",sep="\t",quote=F,row.names=F)




install.packages("VennDiagram")
library(VennDiagram)

TCGA<-read.table("TCGA.txt")
TCGA<-as.matrix(TCGA)
fe<-read.table("Fe.txt")
fe<-as.matrix(fe)

setwd("C:/Users/Administrator/Desktop")
venn.diagram(
  x = list(
    'FerrGenes' =fe ,
    'TCGA_Genes' = TCGA),
  filename = 'intersectGene.png',
  col = "black",
  fill = c("dodgerblue", "goldenrod1"),
  alpha = 0.5,
  cex = 0.8,
  cat.col = 'black',
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.05,
  main = "交集基因数目",
  main.cex = 1.2
)

###################################################
library(VennDiagram)

TCGA<-read.table("Tumor-vs-Normal-diffSig_logFC1_P0.05.xls",sep="\t",header=T)
TCGA<-as.matrix(TCGA$gene)
fe<-read.table("C1-vs-C2-diffSig_logFC1_P0.05.xls",sep="\t",header=T)
fe<-as.matrix(fe$gene)


venn.diagram(
  x = list(
    'DEGs of Tumor-vs-Normal' =TCGA ,
    'DEGs of C1-vs-C2' =fe ),
  filename = 'intersectGene.tif',
 # col = "black",
  fill = c("red", "Blue"),
  #alpha = 0.5,
  resolution = 300, imagetype = "tiff",
  cex = 1,
 # cat.col = 'black',
  #cat.cex = 0.8,
  #cat.fontface = 4,    #main.fontface: 字体样式，比如斜体，粗体等
  margin = 0.05,
  main = "交集基因数目",
  main.cex = 1.2,main.fontface = 2,
  #main.fontfamily = 3,
  #resolution=300,                        #resolution: 输出图形的清晰度，DPI数值
)

A <-intersect(TCGA,fe)
write.table(A,"intersectGenes.txt",quote = F,row.names = F)






