


BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)

######## ######## 使用ALL示例数据#
library(ALL)
data(ALL)
d=exprs(ALL)
d[1:5,1:5]
dim(d)

########## 导入我的数据  ################
dd <- read.table("uniCox_input.txt",sep="\t",header=T,check.names = F,row.names = 1)

unicox <- read.table("uniCox.txt",sep="\t",header=T,check.names = F,row.names = 1)
immune  <- unicox[unicox$pvalue<0.05,]

colnames(dd)
ddd <- dd[,rownames(immune)]     ###提取p<0.05的免疫细胞


dd_t <- as.data.frame(t(ddd))
dd_t[1:5,1:5]


##########  过滤免疫细胞
#data2 <- dd_t[apply(dd_t,1,function(x){sum(is.na(x)) < ncol(dd_t)/2}),]  ##过滤50%缺失值

#c <- rownames(dd_t)[rowSums(dd_t>0) > ncol(dd_t)/2]
#dtt <- dd_t[c,]

##############  不过滤免疫细胞的话就直接往后做

d <- as.matrix(dd_t)


#筛选前5000标准差的基因
#d<- dd_t
#d <- as.matrix(d)
#mads=apply(d,1,mad)  
#d=d[rev(order(mads))[1:nrow(dd_t)],]

#sweep函数减去中位数进行标准化
dt = sweep(d,1, apply(d,1,median,na.rm=T))


#一步完成聚类
library(ConsensusClusterPlus)
#title=tempdir()
dir.create("result_pdf")  #新建文件夹
title="E:\\小逗的工作内容\\项目\\01-大项目\\2021-项目\\24_YQ139-5_BRCA_2021.7.31\\02-分析过程-YQ139-5\\03-ssGSEA-uniCox_一致性聚类\\03-一致性聚类\\result_pdf"
#title=getwd()   #当前路径
results = ConsensusClusterPlus(dt,maxK=10,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",#writeTable=TRUE,
                               seed=1262118388.71279,plot="pdf",tmyPal=NULL)

###d中列-样本，行-features，可以是基因表达矩阵
#maxK：聚类结果中分类的最大数目，必须是整数
#reps:重抽样的次数
#pltem：样品的抽样比例，如pltem=0.8表示采用重抽样方案对样本的80%抽样，经过多次采样，找到稳定可靠的亚组分类
#pFeature：Feature的抽样比例
#clusterAlg：使用的聚类算法，"hc"用于层次聚类
#title：设置生成的文件路径
#distance：计算距离的方法，有pearson，sperson等
#tmyPal：可以指定一致性矩阵使用的颜色，默认使用白-蓝色
#seed：设置随机种子
#plot：不设置时图片结果仅输出屏幕，也可以设置输出为“pdf”,"png"等
#writeTable:若为True,则将一致性矩阵，ICL，log输出到csv文件

#聚类数目K=2，3，4，·····6，采用重抽样方案对样本的80%抽样，经过多次采样，找到稳定可靠的亚组分类。

#输出K=2时的一致性矩阵
results[[2]][["consensusMatrix"]][1:5,1:5]

#hclust选项
results[[2]][["consensusTree"]]   

#样本分类
results[[2]][["consensusClass"]][1:5]  #可以看到各个样品被分到了哪个类别里面去


#计算聚类一致性 (cluster-consensus) 和样品一致性 (item-consensus)
icl <- calcICL(results, title = title, plot = "png")
#writeTable:若为True,则将一致性矩阵，ICL，log输出到csv文件

## 返回了具有两个元素的list，然后分别查看一下
dim(icl[["clusterConsensus"]])
icl[["clusterConsensus"]] 

dim(icl[["itemConsensus"]])

icl[["itemConsensus"]]#[1:5,] 


two <- as.data.frame(results[[3]][["consensusClass"]])
colnames(two) <- "k"

identical(rownames(two),rownames(dd))
k_time <- cbind(dd[,1:2],two)   ###获得了futime，fustst和K


identical(rownames(k_time),rownames(ddd))
k_time_immuneCell <- cbind(k_time,ddd)   ###获得免疫细胞，futime，fustat和k

k_time_immuneCell_1 <- k_time_immuneCell[order(k_time_immuneCell$k),]
k_time_immuneCell_2 <- cbind(id=rownames(k_time_immuneCell_1),k_time_immuneCell_1)

table(k_time_immuneCell_2$k)
write.csv(k_time_immuneCell_2,"样本为3组的结果.csv",row.names=F,quote=F)

##单纯导出样本分组可以用
write.csv(results[[2]][["consensusClass"]],"sample_group2.csv")
write.csv(results[[3]][["consensusClass"]],"sample_group3.csv")


################  亚型间生存分析  #################
k_time_immuneCell_2<- read.csv("样本为3组的结果.csv")

library(survival)
library(survminer)

library(stringr)

#rt=read.table("risk.txt",header=T,sep="\t")
rt <- k_time_immuneCell_2

rt$tile <- c("C")
rt$Subtype <- str_c(rt$tile,rt$k)

table(rt$Subtype)
#write.table(rt,"three_group.txt",sep="\t",quote = F,row.names = F)

rt$futime <- rt$futime/365

diff=survdiff(Surv(futime, fustat) ~Subtype,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ Subtype, data = rt)
surv_median(fit)  #查看高低风险组的中位生存时间
#risk=High_risk 4.068493 
# risk=Low_risk 6.312329

#############  绘制生存曲线

pdf(file="survival_whole.pdf",onefile = FALSE,
    width =7,             #图片的宽度
    height =7)             #图片的高度
ggsurvplot(fit, 
           data=rt,
           conf.int=T,                     #是否展示置信区间
           pval=paste0("p=",pValue),
           pval.size=4,                       #pval字体大小
           risk.table=TRUE,
           legend.labs=c("C1", "C2","C3"),
           legend.title="Subtype",
           xlab="Time(years)",
           break.time.by = 2,            #调整间距
           risk.table.title="",
           #palette=c("red", "blue"),
           palette = c("#E7B800", "#2E9FDF","IndianRed"),
           #palette = c("#E7B800", "#2E9FDF"),
           #palette = c("#F08080", "#00E5EE"),
           risk.table.height=.25,
           risk.table.y.text.col = T, ## risk table左侧是否用对应分层变量的颜色注释
           risk.table.y.text = F, ## 是否展示分层变量的文本，如果不展示，则用颜色条进行展示
           #tables.theme = theme_cleantable(),
           ggtheme = theme_bw() ,# Change ggplot2 theme
           linetype = c('solid', 'solid', 'solid'), surv.median.line = 'hv',)

dev.off()


##################  C1和C2   ####################

library(survival)
library(survminer)

#rt=read.table("risk.txt",header=T,sep="\t")
rt <- k_time_immuneCell_2
rt <- rt[which(rt$k!="3"),]
table(rt$k)

rt$Subtype <- ifelse(rt$k =="1","C1","C2");table(rt$Subtype)

rt$futime <- rt$futime/365

diff=survdiff(Surv(futime, fustat) ~Subtype,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ Subtype, data = rt)
surv_median(fit)  #查看高低风险组的中位生存时间
#risk=High_risk 4.068493 
# risk=Low_risk 6.312329

#############  绘制生存曲线

pdf(file="survival_1和2.pdf",onefile = FALSE,
    width =6,             #图片的宽度
    height =6)             #图片的高度
ggsurvplot(fit, 
           data=rt,
           conf.int=T,                     #是否展示置信区间
           pval=paste0("p=",pValue),
           pval.size=4,                       #pval字体大小
           risk.table=TRUE,
           legend.labs=c("C1", "C2"),
          legend.title="Subtype",
           xlab="Time(years)",
           break.time.by = 2,            #调整间距
           risk.table.title="",
           #palette=c("red", "blue"),
           #palette = c("#E7B800", "#2E9FDF","IndianRed"),
           palette = c("#E7B800", "#2E9FDF"),
           #palette = c("#F08080", "#00E5EE"),
           risk.table.height=.25,
          risk.table.y.text.col = T, ## risk table左侧是否用对应分层变量的颜色注释
          risk.table.y.text = F, ## 是否展示分层变量的文本，如果不展示，则用颜色条进行展示
           #tables.theme = theme_cleantable(),
           ggtheme = theme_bw() ,# Change ggplot2 theme
           linetype = c('solid', 'solid'), surv.median.line = 'hv',)

dev.off()


##################  C1  和 C3   #################
rt <- k_time_immuneCell_2
rt <- rt[which(rt$k!="2"),]
table(rt$k)

rt$Subtype <- ifelse(rt$k =="1","C1","C3");table(rt$Subtype)

rt$futime <- rt$futime/365

diff=survdiff(Surv(futime, fustat) ~Subtype,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ Subtype, data = rt)
surv_median(fit)  #查看高低风险组的中位生存时间
#risk=High_risk 4.068493 
# risk=Low_risk 6.312329

#############  绘制生存曲线

pdf(file="survival_1和3.pdf",onefile = FALSE,
    width =6,             #图片的宽度
    height =6)             #图片的高度
ggsurvplot(fit, 
           data=rt,
           conf.int=T,                     #是否展示置信区间
           pval=paste0("p=",pValue),
           pval.size=4,                       #pval字体大小
           risk.table=TRUE,
           legend.labs=c("C1", "C3"),
           legend.title="Subtype",
           xlab="Time(years)",
           break.time.by = 2,            #调整间距
           risk.table.title="",
           #palette=c("red", "blue"),
           #palette = c("#E7B800", "#2E9FDF","IndianRed"),
           #palette = c("#E7B800", "#2E9FDF"),   #C1和C2
           palette = c("#E7B800", "IndianRed"), #C1和C3
           #palette = c("#F08080", "#00E5EE"),
           risk.table.height=.25,
           risk.table.y.text.col = T, ## risk table左侧是否用对应分层变量的颜色注释
           risk.table.y.text = F, ## 是否展示分层变量的文本，如果不展示，则用颜色条进行展示
           #tables.theme = theme_cleantable(),
           ggtheme = theme_bw() ,# Change ggplot2 theme
           linetype = c('solid', 'solid'), surv.median.line = 'hv',)

dev.off()



##################  C2  和 C3   #################
rt <- k_time_immuneCell_2
rt <- rt[which(rt$k!="1"),]
table(rt$k)

rt$Subtype <- ifelse(rt$k =="2","C2","C3");table(rt$Subtype)

rt$futime <- rt$futime/365

diff=survdiff(Surv(futime, fustat) ~Subtype,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ Subtype, data = rt)
surv_median(fit)  #查看高低风险组的中位生存时间
#risk=High_risk 4.068493 
# risk=Low_risk 6.312329

#############  绘制生存曲线

pdf(file="survival_2和3.pdf",onefile = FALSE,
    width =6,             #图片的宽度
    height =6)             #图片的高度
ggsurvplot(fit, 
           data=rt,
           conf.int=T,                     #是否展示置信区间
           pval=paste0("p=",pValue),
           pval.size=4,                       #pval字体大小
           risk.table=TRUE,
           legend.labs=c("C2", "C3"),
           legend.title="Subtype",
           xlab="Time(years)",
           break.time.by = 2,            #调整间距
           risk.table.title="",
           #palette=c("red", "blue"),
           #palette = c("#E7B800", "#2E9FDF","IndianRed"),
           #palette = c("#E7B800", "#2E9FDF"),   #C1和C2
           #palette = c("#E7B800", "IndianRed"), #C1和C3
           palette = c("#2E9FDF","IndianRed"),   #C2和C3
           #palette = c("#F08080", "#00E5EE"),
           risk.table.height=.25,
           risk.table.y.text.col = T, ## risk table左侧是否用对应分层变量的颜色注释
           risk.table.y.text = F, ## 是否展示分层变量的文本，如果不展示，则用颜色条进行展示
           #tables.theme = theme_cleantable(),
           ggtheme = theme_bw() ,# Change ggplot2 theme
           linetype = c('solid', 'solid'), surv.median.line = 'hv',)

dev.off()
