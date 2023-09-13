


#####   01.提取Tumor中28种免疫基因集的表达矩阵  #############
DEimmune =read.table("ssgseaOut.txt",header=T,sep="\t",check.names=F,row.names = 1)
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

##########  02.导入亚型分组，并和免疫细胞含量merge   ###############

load("clinical_subtype.Rdata")
colnames(merge)

merge1 <- merge[,-c(2:10)]

df<- merge(DEimmune2,merge1,by="id",all=F)

#df <- a[order(a$Subtype),]

###########  03,导入评分文件，并和前边得到的文件合并  ####################
library(dplyr)
score=read.table("scores.txt",sep="\t",check.names=F,row.names=1,header=T)

score$id <- substr(rownames(score),1,12)
score1 <- score %>% select(id, everything())

scores2 <- score1[df$id,]

identical(scores2$id,df$id)
hebing <- merge(df,scores2,by="id",all=F)
hebing <- hebing[order(hebing$Subtype),]

rownames(hebing) <- hebing[,1]
hebing1 <- hebing[,-1]
save(hebing1,file="score_immuneCell_Cluster_TMN.Rdata")
write.csv(hebing1,"score_immuneCell_Cluster_TMN.csv",quote=F)

############### 得到免疫评分，免疫细胞，亚型和临床特征信息的综合文件  ###########

##############      04.绘制多分组热图   #########################
load("score_immuneCell_Cluster_TMN.Rdata")

colnames(hebing1)
gene_matrix <- hebing1[,c(1:28)]
moduleGenes_Exp <- as.data.frame(t(gene_matrix))   ###28种免疫基因集作为表达数据

cluster <- hebing1[,c(29,36:ncol(hebing1))]   #Subtype和ImmuneScore等评分当作分组信息
colnames(cluster)

#cluster$age <- ifelse(cluster$age>"50",">50","<=50")
#cluster$stage <- gsub(' ','',cluster$stage)
#cluster$race <- gsub(' ','_',cluster$race)
colnames(cluster)

#绘制热图
library(pheatmap)
pdf("multiGroup_heatmap.pdf",height=6,width=9)
pheatmap(moduleGenes_Exp, 
         annotation=cluster, 
         #annotation_legend = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3,
         #cutree_col = 3,
         )
dev.off()

################### 05，绘制不同亚型中28种免疫基因集的箱线图  #############################
install.packages("tidyr")
install.packages("ggplot2")
install.packages("tidyverse")
library(tidyverse)
library(tidyr)
library(ggplot2)
require(ggpubr)
require(ggsci)
require(cowplot)
library(RColorBrewer)

#load("score_immuneCell_Cluster_TMN.Rdata")
multi =read.csv("score_immuneCell_Cluster_TMN.csv",sep=",",check.names=F,row.names=1,header=T)
colnames(multi)

hlaexp <- multi[,c(29,1:28)]

boxdata <- gather(hlaexp,HLA,exp,2:ncol(hlaexp))
#boxdata$HLA=gsub("\\.","\\-",boxdata$HLA)
#gather之后的boxdata包含三列信息，表达量信息列，分组信息列，基因名称信息列
boxdata[,'exp'] <- sapply(boxdata[,'exp'],function(x){log2(x+1)})

#group分组列变成factor,boxplot按照L,H排列
boxdata$Subtype <- factor(boxdata$Subtype,levels = c('C1','C2','C3'))

#boxdata <- as.data.frame(apply(boxdata$group,2,function(x) paste('x','risk',sep='-')))
# png('HLA.boxplot.png',width = 650,height = 350)

#display.brewer.all()  #显示所有可用颜色
mycol1 <- brewer.pal(10,'PuOr')
mycol2 <- brewer.pal(10,'BrBG')
mycol3 <- brewer.pal(10,'RdYlBu')

ann_colors =c(mycol3[9],mycol2[9],mycol1[9])
#mycol <- brewer.pal(3,'Dark2')

pdf('Immune_清新.pdf',width =8,height =5)
ggboxplot(boxdata,x='HLA',y='exp',color ='Subtype',
          ylab = 'Relative Expression',
          xlab ='',palette = mycol,  #nejm #jco #jama#npg#d3#aaas
          add ='jitter',add.params = list(size=0.4,shape=18),   #add ='jitter'  shape=16是圆圈
          width=0.6)+
  rotate_x_text(45)+
  #coord_flip()+   #颠倒X轴和Y轴
  theme(axis.text=element_text(size=8),
        axis.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  stat_compare_means(size = 3,aes(group=Subtype),
                     label = 'p.signif',method = 'kruskal.test',hide.ns = F)  ##检验方法可以改
dev.off()

##############     06.绘制各亚型的免疫评分，基质评分和肿瘤纯度的差异   #########################

library(ggplot2)
library(ggpubr)
library(magrittr)
library(ggsignif)

rt <- cluster

###igv,jama
#4组
rt <- rt[order(rt$Subtype),]
pdf("cluster_ImmuneScore.box.pdf",height=6,width=6)
ggboxplot(rt, x= 'Subtype', y='ImmuneScore',
          ylab = "ImmuneScore",xlab = "", color = 'Subtype', palette = "aaas", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Subtype))
dev.off()


pdf("cluster_StromalScore.box.pdf",height=6,width=6)
ggboxplot(rt, x= 'Subtype', y='StromalScore',
          ylab = "StromalScore",xlab = "", color = 'Subtype', palette = "aaas", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Subtype))
dev.off()



pdf("cluster_TumorPurity.box.pdf",height=6,width=6)
ggboxplot(rt, x= 'Subtype', y='TumorPurity',
          ylab = "TumorPurity",xlab = "", color = 'Subtype', palette = "aaas", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Subtype))
dev.off()

















