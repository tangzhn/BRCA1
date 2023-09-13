

##########     01,提取3个基因表达量           ###############################                                    #输入文件

#读取输入文件，并对输入文件处理
rt=read.table("uniq.symbol.txt",sep="\t",header=T,check.names=F)

a <- rt[which(rt$id=="CD274"),]
b <- rt[which(rt$id=="PDCD1"),]
c <- rt[which(rt$id=="CTLA4"),]

CYT <- rbind(a,b,c)

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

write.csv(CYT_3,"免疫检查点分子.csv",quote=F,row.names = F)

##########  02.导入亚型分组和计算几何平均值后的CYT merge   ###############

CYT <- read.table("免疫检查点分子.csv",sep=",",header=T)

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
p <- ggplot(df, aes(x = Subtype, y = CTLA4, color = Subtype)) +
  geom_boxplot(outlier.color = NA) + #隐去箱线图上的异常点
  scale_color_manual(values = mycol) + #自定义配色
  stat_compare_means(#paired = T, #whether you want a paired test
    #两组对比，用wilcox.test或t.test
    method = "kruskal.test", 
    #多组对比，用kruskal.test或anova
    #method = "kruskal.test",
    label.y = max(df$CTLA4)*1.1) + #label的位置
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_line(colour = "black")) + #沿坐标轴显示直线
  xlab("Subtype") + ylab("Expression of CTLA4")+
  #ylim(0,5) + #设置Y轴范围
  guides(color=FALSE) #不显示图例
p

ggsave("Subtype_CTLA4.pdf",height=7,width=7)


# 2，绘制Wikinson点图
p + geom_dotplot(binaxis = "y", #沿y轴堆积，并沿着x轴分组
                 binwidth = 0.5, #最大组距
                 dotsize = 1, #点的大小
                 #如果点太多，两组叠在一起，就需要运行下面这行把它们分开
                 #stackgroups = T, binpositions="all",
                 stackdir = "center")  #数量保持一致的中心堆叠方式

# 3,或者散点图
pdf("dot_Subtype_CTLA4.pdf",height=6,w=6)
p + geom_point(aes(group = Subtype),
               alpha=.3, #点太多，设为透明色，就能看到叠加效果
               size = 2, #点的大小
               position="jitter") #分散
dev.off()







