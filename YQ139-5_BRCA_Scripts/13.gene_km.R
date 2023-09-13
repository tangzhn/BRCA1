
# yy   处理GEO数据
gpl <- read.table('GPL570-55999-58812-88770-20685.txt',header = TRUE,sep = '\t',quote = '', comment.char = '#',check.names = F,fill = T)
#是表示后面的东东都是 标注不要读取
colnames(gpl)


library(dplyr)
probe2id <- gpl %>%
  dplyr::select('ID','Gene Symbol') %>% 
  tidyr::separate_rows('Gene Symbol',sep = ' /// ') %>%
  dplyr::rename(id = 'ID',symbol='Gene Symbol') %>%
  dplyr::filter(symbol != '')


##gene_matrix的获取

probeMatrix <- data.table::fread("GSE88770_series_matrix.txt",data.table = F,skip="ID_REF")
ge <- merge(probe2id,probeMatrix,by.x='id',by.y="ID_REF",all.x=F) 
colnames(ge)[1:20]
gee <- ge[,c(2:ncol(ge))]
boxplot(gee[,90:ncol(gee)],las=2)  #看是否需要标准化
#a <- gee[which(gee$symbol=="IGHG4"),]


##去重
library(limma)
library(estimate)
rt=as.matrix(gee)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]   #把样本表达量数据提取出来
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#data1 <- log2(data+1)   ###标准化数据
data1 <- data
data2 <- data1[rowMeans(data1)>0,]
data <- data2
data <- data.frame(data)

##################### 02.提取25个基因表达矩阵   ###############
aim <- read.table("intersectGenes.txt",sep="\t",header=T)
aim_Exp <- data[aim$x,]

aim_Exp <- na.omit(aim_Exp)   ###去除NA
aim_Exp1 <- data.frame(t(aim_Exp))
aim_Exp2 <- cbind(id=rownames(aim_Exp1),aim_Exp1)

####################  03.导入futime和fustat，与25个基因表达矩阵合并   ##################
time  <- read.table("GSE20685_time.txt",header=T,sep="\t")

rt <- merge(time,aim_Exp2,by="id",all=F)
write.table(rt,"GSE20685_geneSurvival_input.txt",sep="\t",row.names = F,quote=F)


#################### 04，批量生存分析   ###############
library(survival)
#如果以月为单位，除以30；以年为单位，除以365
picDir="05_GSE20685_基因的survival_result"                                               #????ͼƬ????Ŀ¼??????
dir.create(picDir)

pFilter=0.05 
#rt$futime=rt$futime/365                                     #????????Ϊ??λ??????30??????Ϊ??λ??????365
outTab_survival=data.frame()

for(gene in colnames(rt[,4:ncol(rt)])){
  a=rt[,gene]<=median(rt[,gene])
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  outTab_survival=rbind(outTab_survival,cbind(gene=gene,pvalue=pValue))
  pValue=round(pValue,3)
  #pValue=format(pValue, scientific = TRUE)
  
  fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
  summary(fit)
  
  if(pValue<pFilter){
    if(pValue<0.001){
      pValue=signif(pValue,4)
      pValue=format(pValue, scientific = TRUE)
    }else{
      pValue=round(pValue,3)
    }
    
    tiff(file=paste(picDir,"\\",gene,".survival.tiff",sep=""),
         width = 14,            #ͼƬ?Ŀ???
         height =14,            #ͼƬ?ĸ߶?
         units ="cm",
         compression="lzw",
         bg="white",
         res=600)
    plot(fit, 
         lwd=2,
         col=c("#E31A1C","#1F78B4"),
         xlab="Time (year)",
         mark.time=T,
         ylab="Survival rate",
         main=paste(gene,"(p=", pValue ,")",sep="") )
    legend("topright", 
           c("High expression","Low expression"), 
           lwd=2, 
           # col=c("red","blue"))
           col=c("#E31A1C","#1F78B4"))
    dev.off()
  }
}
write.table(outTab_survival,file="05_GSE20685_survival.xls",sep="\t",row.names=F,quote=F)






