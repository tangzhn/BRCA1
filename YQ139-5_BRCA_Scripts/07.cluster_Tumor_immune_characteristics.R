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
gmtFile="APM_IFN.gmt"                                           #GMT文件

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
write.table(ssgseaOut,file="APM_IFN-γ_ssgseaOut.txt",sep="\t",quote=F,col.names=F)

############## 02 和亚型匹配   ###########################################
#####   01.提取Tumor中28种免疫基因集的表达矩阵  #############
DEimmune =read.table("APM_IFN-γ_ssgseaOut.txt",header=T,sep="\t",check.names=F,row.names = 1)
DEimmune1 <- data.frame(t(DEimmune))
colnames(DEimmune1)

###   跳出Tumor样本
DEimmune1$type <- ifelse(substr(rownames(DEimmune1),14,15)  =="01","Tumor","Normal")
table(DEimmune1$type)
DEimmune1 <- subset(DEimmune1,DEimmune1$type=="Tumor")

DEimmune2  <- DEimmune1[,-ncol(DEimmune1)]    
colnames(DEimmune2)

DEimmune2$id <- substr(rownames(DEimmune2),1,12)      ###命名和time一致
DEimmune2 <- DEimmune2[,c(ncol(DEimmune2),1:(ncol(DEimmune2)-1))]  ##也可tarexp <- tarexp %>%  select(id,everything())

##########  02.导入亚型分组，并和APM，和IFN-γ含量merge   ###############

load("clinical_subtype.Rdata")
colnames(merge)

merge1 <- merge[,-c(2:10)]

df<- merge(DEimmune2,merge1,by="id",all=F)
colnames(df)

########得到APM" ,IFN" ,"KEYNOTE.012","KEYNOTE.059" "Subtype" 的文件

##########  03.绘制各亚型下的APM和IFN-γ含量的箱线图   ###############

#开始画图
#用ggplot2画图，用ggpubr算p value。
library(ggpubr)
library(RColorBrewer)
display.brewer.all()  #显示所有可用颜色
#mycol <- brewer.pal(6,'Set1')

mycol <- c("darkgreen", "darkorchid3", "orange") #与分组数量一致 

###1.简单版箱线图  ##########
df <- df[order(df$Subtype),]
p <- ggplot(df, aes(x = Subtype, y = APM, color = Subtype)) +
  geom_boxplot(outlier.color = NA) + #隐去箱线图上的异常点
  scale_color_manual(values = mycol) + #自定义配色
  stat_compare_means(#paired = T, #whether you want a paired test
    #两组对比，用wilcox.test或t.test
    method = "kruskal.test", 
    #多组对比，用kruskal.test或anova
    #method = "kruskal.test",
    label.y = max(df$APM)*1.1) + #label的位置
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_line(colour = "black")) + #沿坐标轴显示直线
  xlab("Subtype") + ylab("Expression of APM")+
  #ylim(0,5) + #设置Y轴范围
  guides(color=FALSE) #不显示图例
p

ggsave("Subtype_APM.pdf",height=7,width=7)


# 2，绘制Wikinson点图
p + geom_dotplot(binaxis = "y", #沿y轴堆积，并沿着x轴分组
                 binwidth = 0.5, #最大组距
                 dotsize = 1, #点的大小
                 #如果点太多，两组叠在一起，就需要运行下面这行把它们分开
                 #stackgroups = T, binpositions="all",
                 stackdir = "center")  #数量保持一致的中心堆叠方式

# 3,或者散点图
pdf("dot_Subtype_APM.pdf",height=6,w=6)
p + geom_point(aes(group = Subtype),
               alpha=.3, #点太多，设为透明色，就能看到叠加效果
               size = 2, #点的大小
               position="jitter") #分散
dev.off()










##########     01,提取CYT的两个基因表达量           ###############################                                    #输入文件

#读取输入文件，并对输入文件处理
rt=read.table("uniq.symbol.txt",sep="\t",header=T,check.names=F)

CYT1 <- rt[which(rt$id=="GZMA"),]
CYT2 <- rt[which(rt$id=="PRF1"),]

identical(colnames(CYT1),colnames(CYT2))
CYT <- rbind(CYT1,CYT2)

boxplot(t(CYT[,2:ncol(CYT)]),las=2)  ###查看是否有0

rownames(CYT) <- CYT[,1]
CYT_1 <- CYT[,-1]
CYT_1t <- data.frame(t(CYT_1))

###   跳出Tumor样本
CYT_1t$type <- ifelse(substr(rownames(CYT_1t),14,15)  =="01","Tumor","Normal")
table(CYT_1t$type)
CYT_2 <- subset(CYT_1t,CYT_1t$type=="Tumor")

CYT_2  <- CYT_2[,-ncol(CYT_2)]    
colnames(CYT_2)

CYT_2$id <- substr(rownames(CYT_2),1,12)      ###命名和time一致
CYT_3 <- CYT_2[,c(ncol(CYT_2),1:(ncol(CYT_2)-1))]  ##也可tarexp <- tarexp %>%  select(id,everything())

write.csv(CYT_3,"CYT.csv",quote=F,row.names = F)

##########  02.导入亚型分组和计算几何平均值后的CYT merge   ###############

CYT <- read.table("CYT_几何平均数文件.txt",sep="\t",header=T)


load("clinical_subtype.Rdata")
colnames(merge)
merge1 <- merge[,-c(2:10)]

df<- merge(CYT,merge1,by="id",all=F)
colnames(df)

########得到CYT" "Subtype" 的文件

##########  03.绘制各亚型下的APM和IFN-γ含量的箱线图   ###############

#开始画图
#用ggplot2画图，用ggpubr算p value。
library(ggpubr)
library(RColorBrewer)
display.brewer.all()  #显示所有可用颜色
#mycol <- brewer.pal(6,'Set1')

mycol <- c("darkgreen", "darkorchid3", "orange") #与分组数量一致 

###1.简单版箱线图  ##########
df <- df[order(df$Subtype),]
p <- ggplot(df, aes(x = Subtype, y = CYT, color = Subtype)) +
  geom_boxplot(outlier.color = NA) + #隐去箱线图上的异常点
  scale_color_manual(values = mycol) + #自定义配色
  stat_compare_means(#paired = T, #whether you want a paired test
    #两组对比，用wilcox.test或t.test
    method = "kruskal.test", 
    #多组对比，用kruskal.test或anova
    #method = "kruskal.test",
    label.y = max(df$CYT)*1.1) + #label的位置
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_line(colour = "black")) + #沿坐标轴显示直线
  xlab("Subtype") + ylab("Cytotoxic activity scores(CYT)")+
  #ylim(0,5) + #设置Y轴范围
  guides(color=FALSE) #不显示图例
p

ggsave("Subtype_CYT.pdf",height=7,width=7)


# 2，绘制Wikinson点图
p + geom_dotplot(binaxis = "y", #沿y轴堆积，并沿着x轴分组
                 binwidth = 0.5, #最大组距
                 dotsize = 1, #点的大小
                 #如果点太多，两组叠在一起，就需要运行下面这行把它们分开
                 #stackgroups = T, binpositions="all",
                 stackdir = "center")  #数量保持一致的中心堆叠方式

# 3,或者散点图
pdf("dot_Subtype_CYT.pdf",height=6,w=6)
p + geom_point(aes(group = Subtype),
               alpha=.3, #点太多，设为透明色，就能看到叠加效果
               size = 2, #点的大小
               position="jitter") #分散
dev.off()










##########     01,提取TIS的18个基因表达量           ###############################                                    #输入文件

#读取输入文件，并对输入文件处理
rt=read.table("uniq.symbol.txt",sep="\t",header=T,check.names=F)
rownames(rt) <- rt$id
rt1 <- rt[,-1]

aim=read.table("TIS_geneList.txt",sep="\t",header=T,check.names=F)

a <- intersect(rownames(rt1),aim$id)
df <- rt1[a,]

TIS <- data.frame(t(df))

###   跳出Tumor样本
TIS$type <- ifelse(substr(rownames(TIS),14,15)  =="01","Tumor","Normal")
table(TIS$type)

TIS1 <- subset(TIS,TIS$type=="Tumor")

TIS2  <- TIS1[,-ncol(TIS1)]    
colnames(TIS2)

TIS2$id <- substr(rownames(TIS2),1,12)      ###命名和time一致
TIS3 <- TIS2[,c(ncol(TIS2),1:(ncol(TIS2)-1))]  ##也可tarexp <- tarexp %>%  select(id,everything())
write.csv(TIS3,"15个基因的表达矩阵.csv",quote=F,row.names = F)

#################    算15个基因的平均数  #############
TI=read.csv("15个基因的表达矩阵.csv",sep=",",header=T,check.names=F)

colnames(TI)
TI$TIS_mean <- rowMeans(TI[,2:ncol(TI)])
write.csv(TI,"TIS.csv",quote=F,row.names = F)

boxplot(TI[,2:ncol(TI)],las=2)

TIS111 <- TI[,c(1,ncol(TI))]   ###样本和TIS提取出来

##########  02.导入亚型分组和计算出的TIS  merge   ###############

load("clinical_subtype.Rdata")
colnames(merge)
merge1 <- merge[,-c(2:10)]

df<- merge(TIS111,merge1,by="id",all=F)
colnames(df)

########得到CYT" "Subtype" 的文件

##########  03.绘制各亚型下的APM和IFN-γ含量的箱线图   ###############

#开始画图
#用ggplot2画图，用ggpubr算p value。
library(ggpubr)
library(RColorBrewer)
display.brewer.all()  #显示所有可用颜色
#mycol <- brewer.pal(6,'Set1')

mycol <- c("darkgreen", "darkorchid3", "orange") #与分组数量一致 

###1.简单版箱线图  ##########
df <- df[order(df$Subtype),]
p <- ggplot(df, aes(x = Subtype, y = TIS_mean, color = Subtype)) +
  geom_boxplot(outlier.color = NA) + #隐去箱线图上的异常点
  scale_color_manual(values = mycol) + #自定义配色
  stat_compare_means(#paired = T, #whether you want a paired test
    #两组对比，用wilcox.test或t.test
    method = "kruskal.test", 
    #多组对比，用kruskal.test或anova
    #method = "kruskal.test",
    label.y = max(df$TIS_mean)*1.1) + #label的位置
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_line(colour = "black")) + #沿坐标轴显示直线
  xlab("Subtype") + ylab("Tumor Inflammation Signature(TIS)")+
  #ylim(0,5) + #设置Y轴范围
  guides(color=FALSE) #不显示图例
p

ggsave("Subtype_TIS.pdf",height=7,width=7)


# 2，绘制Wikinson点图
p + geom_dotplot(binaxis = "y", #沿y轴堆积，并沿着x轴分组
                 binwidth = 0.5, #最大组距
                 dotsize = 1, #点的大小
                 #如果点太多，两组叠在一起，就需要运行下面这行把它们分开
                 #stackgroups = T, binpositions="all",
                 stackdir = "center")  #数量保持一致的中心堆叠方式

# 3,或者散点图
pdf("dot_Subtype_TIS.pdf",height=6,w=6)
p + geom_point(aes(group = Subtype),
               alpha=.3, #点太多，设为透明色，就能看到叠加效果
               size = 2, #点的大小
               position="jitter") #分散
dev.off()










##########              ###############################                                    #输入文件
aim=read.table("ssgseaOut.txt",sep="\t",header=T,check.names=F,row.names = 1)
aim1 <- data.frame(t(aim))
#aim2 <- cbind(id=rownames(aim1),aim1)


###   跳出Tumor样本
aim1$type <- ifelse(substr(rownames(aim1),14,15)  =="01","Tumor","Normal")
table(aim1$type)

TIS1 <- subset(aim1,aim1$type=="Tumor")
TIS2  <- TIS1[,-ncol(TIS1)]    
colnames(TIS2)

TIS2$id <- substr(rownames(TIS2),1,12)      ###命名和time一致
TIS3 <- TIS2[,c(ncol(TIS2),1:(ncol(TIS2)-1))]  ##也可tarexp <- tarexp %>%  select(id,everything())
#write.csv(TIS3,".csv",quote=F,row.names = F)

#################    算15个基因的平均数  #############
TI=read.csv("15个基因的表达矩阵.csv",sep=",",header=T,check.names=F)

colnames(TI)
TI$TIS_mean <- rowMeans(TI[,2:ncol(TI)])
write.csv(TI,"TIS.csv",quote=F,row.names = F)

boxplot(TI[,2:ncol(TI)],las=2)

TIS111 <- TI[,c(1,ncol(TI))]   ###样本和TIS提取出来

##########  02.导入亚型分组和计算出的TIS  merge   ###############

load("clinical_subtype.Rdata")
colnames(merge)
merge1 <- merge[,-c(2:10)]

df<- merge(TIS3,merge1,by="id",all=F)
colnames(df)

########得到CYT" "Subtype" 的文件

##########  03.绘制各亚型下的APM和IFN-γ含量的箱线图   ###############

#开始画图
#用ggplot2画图，用ggpubr算p value。
library(ggpubr)
library(RColorBrewer)
display.brewer.all()  #显示所有可用颜色
#mycol <- brewer.pal(6,'Set1')

mycol <- c("darkgreen", "darkorchid3", "orange") #与分组数量一致 

###1.简单版箱线图  ##########
df <- df[order(df$Subtype),]
p <- ggplot(df, aes(x = Subtype, y = TIL, color = Subtype)) +
  geom_boxplot(outlier.color = NA) + #隐去箱线图上的异常点
  scale_color_manual(values = mycol) + #自定义配色
  stat_compare_means(#paired = T, #whether you want a paired test
    #两组对比，用wilcox.test或t.test
    method = "kruskal.test", 
    #多组对比，用kruskal.test或anova
    #method = "kruskal.test",
    label.y = max(df$TIL)*1.1) + #label的位置
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_line(colour = "black")) + #沿坐标轴显示直线
  xlab("Subtype") + ylab("Tumor Infiltrating Lymphocytos(TILs)")+
  #ylim(0,5) + #设置Y轴范围
  guides(color=FALSE) #不显示图例
p

ggsave("Subtype_TILs.pdf",height=7,width=7)


# 2，绘制Wikinson点图
p + geom_dotplot(binaxis = "y", #沿y轴堆积，并沿着x轴分组
                 binwidth = 0.5, #最大组距
                 dotsize = 1, #点的大小
                 #如果点太多，两组叠在一起，就需要运行下面这行把它们分开
                 #stackgroups = T, binpositions="all",
                 stackdir = "center")  #数量保持一致的中心堆叠方式

# 3,或者散点图
pdf("dot_Subtype_TILs.pdf",height=6,w=6)
p + geom_point(aes(group = Subtype),
               alpha=.3, #点太多，设为透明色，就能看到叠加效果
               size = 2, #点的大小
               position="jitter") #分散
dev.off()







