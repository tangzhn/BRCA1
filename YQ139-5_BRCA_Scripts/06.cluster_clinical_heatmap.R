#    箱线图  ########
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

hlaexp<-read.table("three_group.txt",sep="\t",header=T,check.names = F)
colnames(hlaexp)
hlaexp <- hlaexp[,-c(2,3,4,14)]

time <- read.table("clinical.txt",sep="\t",header=T,check.names = F)
time <- time[,-2:-3]


merge <- merge(hlaexp,time,by="id",all=F)
colnames(merge)

save(merge,file="clinical_subtype.Rdata")


load("clinical_subtype.Rdata")
rownames(merge) <- merge$id
merge <- merge[,-1]


df <- merge[,c(8:ncol(merge))]
df <- df[order(df$Subtype),]
colnames(df)

df2 <- as.data.frame(t(df[,1:2]))

cluster <- df[,c(4:ncol(df),3)]
cluster <- cluster[order(cluster$Subtype),]
cluster$age <- ifelse(cluster$age>"50","dayu50","xiaoyudengyu50")
cluster$stage <- gsub(' ','',cluster$stage)
cluster$race <- gsub(' ','_',cluster$race)
colnames(cluster)

#绘制热图
library(pheatmap)
#pdf("multiGroup_heatmap--3.pdf",height=8,width=8)
pheatmap(df2, 
         annotation=cluster, 
         #annotation_legend = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3,
         cutree_col = 3)
#dev.off()


annotation_col = data.frame(Age = factor(cluster$age),
                            M = factor(cluster$M),
                            N = factor(cluster$N),
                            T = factor(cluster$T),
                            Stages = factor(cluster$stage),
                            Race = factor(cluster$race),
                            Subtype = factor(cluster$Subtype))  ##分组，legend名称：Type,这里可以调整显示顺序
rownames(annotation_col) = rownames(cluster)
#annotation_col <- annotation_col[,c(2,3,7)]

unique(cluster$M);unique(cluster$N);unique(cluster$stage);unique(cluster$T);unique(cluster$race)
library(RColorBrewer)
#display.brewer.all()  #显示所有可用颜色

#给分组设置颜色.注意格式
ann_age <-  brewer.pal(4,"PuRd")
#ann_age <- c('#FBB4A1','#B3CDE3')
#mycol1 <- brewer.pal(6,'Set1')
#ann_M <- c('#7FC97F','#BEAED4')
ann_M <- brewer.pal(6,'PuBuGn')
ann_N <- brewer.pal(4,'YlGnBu')
ann_stage <- brewer.pal(10,'BrBG')
ann_T <- brewer.pal(5,'YlOrRd')
ann_race <- brewer.pal(5,'Dark2')
ann_Subtype <- brewer.pal(6,'Set2')

ann_colors = list(Age = c(dayu50=ann_age[3],xiaoyudengyu50=ann_age[1]),
                  M = c(M0=ann_M[5],M1=ann_M[2]),
                  N = c(N0=ann_N[1],N1=ann_N[2],N2=ann_N[3],N3=ann_N[4]),
                  T = c(T1=ann_T[1],T2=ann_T[3],T3=ann_T[4],T4=ann_T[6]),
                  Stages = c(StageI=ann_stage[7],StageII=ann_stage[8],StageIII=ann_stage[9],StageIV=ann_stage[10]),
                  Race = c(black_or_african_american=ann_race[1],white=ann_stage[5],asian=ann_T[2],american_indian_or_alaska_native=ann_race[2],not_reported=ann_race[5]),
                  Subtype = c(C1=ann_Subtype[2],C2=ann_Subtype[3],C3=ann_Subtype[5]))

library(pheatmap)
pdf("multiGroup_heatmap--3.pdf",height=8,width=8)
identical(rownames(annotation_col),colnames(df2))


pheatmap(df2, 
         annotation=annotation_col,
         annotation_colors = ann_colors,##
         #annotation_legend = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3,
         main = "C1(n=424) C2(n=500) C3(n=115)",
         cutree_cols = 3,
         )
dev.off()





rownames(hlaexp) <- hlaexp[,1]
#rownames(hlaexp) <- hlaexp$ID 
hlaexp=hlaexp[,-1]












boxdata <- gather(hlaexp,HLA,exp,2:ncol(hlaexp))
#boxdata$HLA=gsub("\\.","\\-",boxdata$HLA)
#gather之后的boxdata包含三列信息，表达量信息列，分组信息列，基因名称信息列
boxdata[,'exp'] <- sapply(boxdata[,'exp'],function(x){log2(x+1)})

#group分组列变成factor,boxplot按照L,H排列
boxdata$group <- factor(boxdata$group,levels = c('Control','Degenerated'))

#boxdata <- as.data.frame(apply(boxdata$group,2,function(x) paste('x','risk',sep='-')))


# png('HLA.boxplot.png',width = 650,height = 350)

pdf('4-gene.pdf',width =7,height = 5)
ggboxplot(boxdata,x='HLA',y='exp',color ='group',
          ylab = 'Relative Expression',
          xlab ='',palette = "d3",  #nejm #jco #jama#npg#d3#aaas
          add ='jitter',add.params = list(size=1,shape=16),   #add ='jitter'
          width=0.6)+
  #rotate_x_text(45)+
  #coord_flip()+   #颠倒X轴和Y轴
  theme(axis.text=element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))+
  stat_compare_means(size = 4,aes(group=group),
                     label = 'p.signif',method = 'wilcox.test',hide.ns = F)  ##检验方法可以改
dev.off()
