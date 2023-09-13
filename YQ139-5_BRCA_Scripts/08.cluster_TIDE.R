if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")  #, version = "3.8"

library("limma")

#setwd("C:\\Users\\Administrator\\Desktop\\11_risk_treatment")
rt<-read.table("uniq.symbol.txt",sep="\t",header=T,row.names=1,check.names = F)
rt[1:3,1:3]

rt<-as.matrix(rt)
rt<-log2(rt+1)

rt1<-mean(rt)
rt<-rt-rt1

#rt$HLA=gsub("\\.","\\-",boxdata$HLA)
write.table(rt,"matrix_normlize.txt",sep="\t",quote=F,row.names=T)




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

#################    01，导入TIDE结果  #############
TIS=read.csv("TIDE_result.csv",sep=",",header=T,check.names=F,row.names = 1)

###   跳出Tumor样本
TIS$type <- ifelse(substr(rownames(TIS),14,15)  =="01","Tumor","Normal")
table(TIS$type)

TIS1 <- subset(TIS,TIS$type=="Tumor")
TIS2  <- TIS1[,-ncol(TIS1)]    
colnames(TIS2)

TIS2$id <- substr(rownames(TIS2),1,12)      ###命名和time一致
TIS3 <- TIS2[,c(ncol(TIS2),1:(ncol(TIS2)-1))]  ##也可tarexp <- tarexp %>%  select(id,everything())
colnames(TIS3)

TIDE <- TIS3[,c(1,4)]

##########  02.导入亚型分组和计算出的TIS  merge   ###############

load("clinical_subtype.Rdata")
colnames(merge)
merge1 <- merge[,-c(2:10)]

df<- merge(TIDE,merge1,by="id",all=F)
colnames(df)

########得到TIDE "Subtype" 的文件

##########  03.绘制各亚型下的TIDE的箱线图   ###############

#开始画图
#用ggplot2画图，用ggpubr算p value。
library(ggpubr)
library(RColorBrewer)
display.brewer.all()  #显示所有可用颜色
#mycol <- brewer.pal(6,'Set1')

mycol <- c("darkgreen", "darkorchid3", "orange") #与分组数量一致 

###1.简单版箱线图  ##########
df <- df[order(df$Subtype),]
p <- ggplot(df, aes(x = Subtype, y = TIDE, color = Subtype)) +
  geom_boxplot(outlier.color = NA) + #隐去箱线图上的异常点
  scale_color_manual(values = mycol) + #自定义配色
  stat_compare_means(#paired = T, #whether you want a paired test
    #两组对比，用wilcox.test或t.test
    method = "kruskal.test", 
    #多组对比，用kruskal.test或anova
    #method = "kruskal.test",
    label.y = max(df$TIDE)*1.5) + #label的位置
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_line(colour = "black")) + #沿坐标轴显示直线
  xlab("Subtype") + ylab("TIDE")+
  #ylim(0,5) + #设置Y轴范围
  guides(color=FALSE) #不显示图例
p

ggsave("Subtype_TIDE.pdf",height=7,width=7)


# 2，绘制Wikinson点图
p + geom_dotplot(binaxis = "y", #沿y轴堆积，并沿着x轴分组
                 binwidth = 0.5, #最大组距
                 dotsize = 1, #点的大小
                 #如果点太多，两组叠在一起，就需要运行下面这行把它们分开
                 #stackgroups = T, binpositions="all",
                 stackdir = "center")  #数量保持一致的中心堆叠方式

# 3,或者散点图
pdf("dot_Subtype_TIDE.pdf",height=6,w=6)
p + geom_point(aes(group = Subtype),
               alpha=.3, #点太多，设为透明色，就能看到叠加效果
               size = 2, #点的大小
               position="jitter") #分散
dev.off()

########################  violin   #########################

display.brewer.all()  #显示所有可用颜色

mycol1 <- brewer.pal(6,'YlGnBu')
mycol2 <- brewer.pal(6,'Reds')
mycol3 <- brewer.pal(6,'BuPu')

mycol <- c(mycol1[4],mycol2[4],mycol3[4]) #与分组数量一致 
mycol4444 <- brewer.pal(6,'Dark2') #与分组数量一致 

p <- ggplot(df, aes(x = Subtype, y = TIDE, fill = Subtype)) +
  geom_violin(outlier.color = NA) + #隐去箱线图上的异常点
  #scale_color_manual(values = mycol) + #自定义配色
  scale_fill_manual(values = mycol4444) + #自定义配色
  stat_compare_means(#paired = T, #whether you want a paired test
    #两组对比，用wilcox.test或t.test
    #method = "kruskal.test", 
    #多组对比，用kruskal.test或anova
    #method = "kruskal.test",
    label.y = max(df$TIDE)*1.5) + #label的位置
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_line(colour = "black")) + #沿坐标轴显示直线
  xlab("Subtype") + ylab("TIDE")+
  #ylim(0,5) + #设置Y轴范围
  guides(color=FALSE) #不显示图例
p
ggsave("Subtype_TIDE_violin_实心.pdf",height=7,width=7)



