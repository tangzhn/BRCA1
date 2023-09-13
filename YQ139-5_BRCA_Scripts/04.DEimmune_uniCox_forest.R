

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
colnames(DEimmune2)

###  筛选出差异免疫细胞  ########
DEimmune3 <- DEimmune2[,-c(5,9,16,22,26)]  ##筛选差异免疫细胞


#############  导入time文件，与免疫细胞数据merge
time=read.table("time.txt",header=T,sep="\t",check.names=F,)

rtt <- merge(time,DEimmune3,by = 'id')
#colnames(rtt) <- gsub(" ",".",colnames(rtt))
write.table(rtt,"uniCox_input.txt",sep="\t",row.names = F,quote=F)




#######  剔除50%以上样本不表达的免疫极阴极，不剔除就不运行  #######

#cell <- colnames(rtt)[4:ncol(rtt)][sapply(4:ncol(rtt), function(i){
#  quantile(rtt[,i],0.5) > 0} )]
#cell

rt <- rtt[,c("id","futime","fustat", cell)]

c <- colnames(rtt)[4:ncol(rtt)][colSums(rtt[,4:ncol(rtt)]>0)>931*5/10]
data.frame(c)
rt <- rtt[,c("id","futime","fustat",c)]


############    单因素Cox回归分析   ###################
library(survival)

rt <- rtt
outTab=data.frame()
for(i in colnames(rt[,4:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     coef=coxSummary$coefficients[,"coef"],  
                     HR=coxSummary$conf.int[,"exp(coef)"],  #View(coxSummary$conf.int)
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)



library(survminer)
names(rt)[names(rt) == 'T.cells.regulatory.(Tregs)'] <- 'T.cells.regulatory'



############绘制森林图函数----一般预后分析用这个   ############

bioForest=function(coxFile=null,color=null,forestFile=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")   #本来的代码
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width = 9,height = 7)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.7
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'P value',cex=text.cex,font=4,adj=1)  #font设置文字格式，1默认，2，加粗，3斜体，4加粗+斜体
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'HR(95% CI)',cex=text.cex,font=2,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(min(as.numeric(hrLow)),max(as.numeric(hrLow),as.numeric(hrHigh)))  ###改HR的起始位置坐标
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="black",lwd=1)  ###横着那根线
  abline(v=1,col="black",lty=2,lwd=1)
 # boxcolor = ifelse(as.numeric(hr) > 1, color, color)
  boxcolor = ifelse(as.numeric(hr) > 1, "#DC143C", "DodgerBlue")
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)  #pch=15是方块
  axis(1)
  dev.off()
}

bioForest(coxFile="uniCox.txt",color="DodgerBlue",forestFile="uniForest.pdf")


####颜色DarkTurquoise





##############   第二种森林图   ################

covariates <- colnames(rt[,4:ncol(rt)])

colnames(rt[4:ncol(rt)])
#names(rt)[names(rt) == 'T.cells.regulatory.(Tregs)'] <- 'T.cells.regulatory'

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(futime, fustat)~', x)))
#对每一个特征做cox回归分析
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = rt)})

#提取HR，95%置信区间和p值
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         #获取HR
                         HR <-signif(x$coef[2], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

#转换成数据框，并转置
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <-as.data.frame(res,stringsAsFactors=F)

#对HR (95% CI for HR)做处理，得到HR和low .95和high .95
#当然也可以改计算univ_results这一步的代码，不要将HR和CI贴起来
############################################################
HR=gsub("[\\(\\)]","",res$`HR (95% CI for HR)`)
HR=gsub("-"," ",HR)
HR=as.data.frame(do.call(cbind,strsplit(HR," ")),stringsAsFactors=F)
names(HR)=rownames(res)


#################################
#开始绘图，直接保存到pdf文件中
#################################
pdf(file="univariate_forest.pdf",width=7)
#左边和右边边距稍微留多一点来写变量名称，pvalue和HR
par(mar=c(5,6,4,13))
#先用小方块画出HR
plot(as.numeric(HR[1,]),1:dim(HR)[2],
     pch=15,cex=2,col="DodgerBlue",bty='n',yaxt='n',ylab=NA,xlab="Hazard Ratio",
     xlim=range(as.numeric(unlist(HR)))
)
#添加中线
abline(v=1,col="grey",lwd=2,lty=1)    #lty=2 是虚线

for(i in 1:ncol(HR)){
  x=as.numeric(HR[2:3,i])
  #循环画出CI
  lines(x,c(i,i),col="DodgerBlue",lwd=3)  #lwd对线加粗
  #添加变量名
  text(0.2,i,rownames(res)[i],xpd=T,adj = c(0,0))
  #添加p值
  text(2.1,i,as.numeric(res[i,1]),xpd=T,adj = c(0,0))
  #添加HR和CI
  text(2.7,i,as.character(res[i,2]),xpd=T,adj = c(0,0))
}
#添加标题
text(2.1,ncol(HR)+0.5,"P value",xpd=T,adj = c(0,0),font=3)
text(2.7,ncol(HR)+0.5,"HR(95% CI)",xpd=T,adj = c(0,0),font=1)
dev.off()


##########  差异免疫细胞KM生存分析    ##################################
library(survival)
rt=read.table("uniCox_input.txt",header=T,sep="\t",check.names=F)      #读取文件
#如果以月为单位，除以30；以年为单位，除以365
picDir="survival_result"                                               #????ͼƬ????Ŀ¼??????
dir.create(picDir)

pFilter=0.05                                                #????ͼ??̫?ֻ࣬??pС??0.05?Ļ?????ͼ
rt$futime=rt$futime/365                                     #????????Ϊ??λ??????30??????Ϊ??λ??????365
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
         col=c("#1F78B4","#E31A1C"),
         xlab="Time (year)",
         mark.time=T,
         ylab="Survival rate",
         main=paste(gene,"(p=", pValue ,")",sep="") )
    legend("topright", 
           c("High expression","Low expression"), 
           lwd=2, 
          # col=c("red","blue"))
           col=c("#1F78B4","#E31A1C"))
    dev.off()
  }
}
write.table(outTab_survival,file="survival.xls",sep="\t",row.names=F,quote=F)









