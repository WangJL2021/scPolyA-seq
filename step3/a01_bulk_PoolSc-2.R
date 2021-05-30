library(ggplot2)
library(Cairo)

setwd("/home/wangjl/data/apa/191111Figure/f1/bulk")
getwd()

##############
# load data

#(1) load bulk
bulk1=read.table("/home/wangjl/data/apa/bulk/htseq/bulkP1_R2.bam.count", row.names = 1,header = F)
colnames(bulk1)="bulkP1"
dim(bulk1) #57792     1
head(bulk1)
#
bulk2=read.table("/home/wangjl/data/apa/bulk/htseq/bulkP2_R2.bam.count", row.names = 1,header = F)
colnames(bulk2)="bulkP2"
dim(bulk2)  #57792     1
head(bulk2)
#
bulkAll0=cbind(bulk1,bulk2)
head(bulkAll0)
dim(bulkAll0) #57792     2

# (2) load sc
sc0=read.csv("/home/wangjl/data/apa/190530Mix/BC_HeLa.225cells.count.V3.csv", row.names = 1)
head(sc0[,1:4]) 
dim(sc0) #18679   225




#########################
head(bulk1)
dim(bulk1)
dim(sc0)

sc1=apply(sc0, 1, sum)
head(sc1)

g1=intersect(rownames(bulk1), names(sc1))
length(g1)


#1. add at raw counts
p1.1=(function(){
  bk1=bulk1[g1,]
  head(bk1)
  
  sc1=sc1[g1]
  head(sc1)
  
  # add counts, no scale
  p1=cor(bk1, unname(sc1))
  p2=cor(log2(1+bk1), log2(1+unname(sc1)) )
  c(p1, p2) #0.678
})()
p1.1 #0.678 0.913246





#2. add at normalized counts.

## method1: if subset first, then normalize?
p2.1=(function(){
  b1=bulk1[g1,]
  s1=sc0[g1,]
  #
  b1.norm=b1/sum(b1)*1e6
  s1.norm=apply(s1, 2, function(x){
    x/sum(x)*1e6
  })
  s1.norm.1=apply(s1.norm, 1, sum)/255
  
  p1=cor(b1.norm, s1.norm.1)
  p2=cor(log2(1+b1.norm), log2(1+s1.norm.1) )
  plot(log2(1+b1.norm),log2(1+s1.norm.1) )
  c(p1, p2)
})()
p2.1 #0.667, 0.954


## method2: if normalize, then subset
p2.2=(function(){
  b1=bulk1
  s1=sc0
  #
  b1.norm=b1/sum(b1)*1e6
  s1.norm=apply(s1, 2, function(x){
    x/sum(x)*1e6
  })
  s1.norm.1=apply(s1.norm, 1, sum)/255
  
  # subset
  b2=b1.norm[g1,"bulkP1"]
  s2=s1.norm.1[g1]
  
  p1=cor(b2, s2)
  p2=cor(log2(1+b2), log2(1+s2) )
  plot(log2(1+b2),log2(1+s2) )
  
  library(ggplot2)
  g=ggplot(NULL, aes(log2(1+b2), log2(1+s2) ) )+theme_classic()+
    geom_point(size=0.02, alpha=0.1)+
    geom_text(aes(x=0,y=14,label=paste0("r=", round(p2,3) ) ),hjust=0 )+  
    labs(x="Tag of bulk data (log2)", y="Tag of pooled single cells (log2)",
         title="Gene expresion in log2(cpm+1)")
  CairoPDF('bulk_PoolSC.cor-3.pdf', width=3.3,height=3.3)
  print(g)
  dev.off()
  
  c(p1, p2)
})()
p2.2 #0.667, 0.957




#3. add at normal, log2 counts
p3.1=(function(){
  b1=bulk1
  s1=sc0
  #
  b1.norm=b1/sum(b1)*1e6
  s1.norm=apply(s1, 2, function(x){
    x/sum(x)*1e6
  })
  s1.norm.1=apply(log2(1+s1.norm), 1, sum)/255
  # subset
  b2=b1.norm[g1,"bulkP1"]
  s2=s1.norm.1[g1]
  
  p1=cor(b2, s2)
  p2=cor(log2(1+b2), log2(1+s2) )
  plot(log2(1+b2),log2(1+s2) )
  
  c(p1, p2)
})()
p3.1 #0.368, 0.9446
