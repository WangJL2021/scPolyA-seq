#c01_check_sigCor_expressionLevel.R

setwd("/data/jinwf/wangjl/apa/20200701Fig/f6/result/exp_BC/")
getwd()


###########
# load data
#(1) cell info
cellInfo=read.table("/data/jinwf/wangjl/apa/20200701Fig/f3/cell_cycle/cellInfo.V6.txt")
dim(cellInfo) #222   8
head(cellInfo)
table(cellInfo$cellType)

cid.BC=rownames(cellInfo[which( substr(cellInfo$cellType,1,2)=='BC'),]); length(cid.BC) #165

#(2) add RNA read counts for each genes
rnaM=read.csv('/data/jinwf/wangjl/apa/20200701Fig/f2/BC_HeLa/BC_HeLa.222cells.count.V4.csv', row.names = 1)
dim(rnaM) #18662 222
rnaM[1:10,1:10]
#
rnaM.logcpm=as.data.frame(apply(rnaM, 2, function(x){
  log2(x/sum(x)*1e6+1)
}))
dim(rnaM.logcpm)
rnaM.logcpm[1:10,1:10]

rnaM.logcpm.BC=rnaM.logcpm[, cid.BC]
dim(rnaM.logcpm.BC) #18662   165




###########
# load gene set
neg=readLines("../cor_BC/MDA-MB-468_02-Spearman_cor_gDPAU_RNA.neg_high.gene.txt")
length(neg) #595

pos=readLines("../cor_BC/MDA-MB-468_02-Spearman_cor_gDPAU_RNA.pos_high.gene.txt")
length(pos) #222

meanExp=as.data.frame(apply(rnaM.logcpm.BC[c(neg, pos),], 1, mean))
colnames(meanExp)="meanLogExp"
head(meanExp)
meanExp$group="pos"
meanExp[neg,'group']="neg"
table(meanExp$group)

# save
write.table(meanExp, "meanExp_sigCor_neg_vs_pos_BC.df.txt")




#############
# plot
# transfer p value for ggplot2 plotting
transPvalue2=function(p0){
  #p0=1.2e-52; p0
  if(p0<0.01){
    p0=formatC(p0, format="e", digits=2);#p0 #"1.2e-52"
    # split by e
    p0=strsplit(p0, "e")[[1]]; #p0 #[1] "1.2" "-52"
    # plot
    label = paste0( 'italic(P)~"="~',p0[1],'~"\u00d7"~10^',p0[2])
  }else{
    p0=round(p0,3); #p0 #保留2位
    label = paste0( 'italic(P)~"=" ~', p0 )
  }
  return(label)
}

p1=t.test(meanExp[which(meanExp$group=="pos"),]$meanLogExp, 
          meanExp[which(meanExp$group=="neg"),]$meanLogExp)
p1$p.value



pdf("expLevel_sigCor_genes.pdf", width=2.4, height=3.8)
boxplot(meanExp$meanLogExp ~ meanExp$group,
        xaxt='n', ylim=c(0,14.5), las=2, 
        cex=0.5, cex.main=0.7,
        col=c('blue', 'red'),
        xlab="", ylab="Average expression(log2)",
        main="MDA-MB-468 cell line\nGroup by APA and RNA correlation")
text(1.5, 14, labels = parse(text=transPvalue2(p1$p.value)), cex=0.8 )

text(x=c(1, 2), y=-1, labels=c("Negative", "Positive"),
     cex=1, xpd=TRUE, adj=1, srt=60)
dev.off()
# end 




table(meanExp$group)
#neg pos 
#595 222 



#####
library(ggplot2)
ggplot(meanExp, aes(group, meanLogExp))+
  geom_boxplot()+
  #geom_jitter()+
  theme_bw()

