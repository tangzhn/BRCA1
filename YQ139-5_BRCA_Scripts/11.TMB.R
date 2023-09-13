if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("maftools")


library(maftools)
#setwd("C:\\Users\\lexb4\\Desktop\\tcgaTMB\\03.maftools")
maf = read.maf(maf = 'C1_out.maf')
maf
write.mafSummary(maf = maf, basename = 'C1-maf')


pdf(file="C3_summary.pdf",width=7,height=6)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf(file="C3_waterfall.pdf",width=7,height=6)
oncoplot(maf = maf, top = 20, fontSize = 0.6 ,showTumorSampleBarcodes = F )
dev.off()

pdf(file="C3_interaction.pdf",width=7,height=6)
somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.01))
dev.off()

###可以使用 oncostrip 函数展示特定基因在样本中的突变情况，
###此处查看肝癌中关注较多的'TP53','CTNNB1', 'ARID1A'三个基因，如下：
oncostrip(maf = maf, genes = c('TP53','CTNNB1', 'ARID1A'))



pdf(file="C1_Genecloud.pdf",width=7,height=6)
geneCloud(input = maf, minMut = 5)
dev.off()



######??????ѧ??: http://study.163.com/u/biowolf
######??????ѧ??: https://shop119322454.taobao.com
######??????ѧ??: http://www.biowolf.cn/
######???????䣺2740881706@qq.com
######????΢??: seqBio
######QQȺ:  259208034
