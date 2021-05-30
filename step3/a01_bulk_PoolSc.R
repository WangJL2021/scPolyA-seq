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



##############
# gene order
genes=intersect(rownames(bulkAll0), rownames(sc0))
length(genes) # 18679

countSc=sc0[genes,]
countBulk=bulkAll0[genes,]
dim(countSc) #18679   225
dim(countBulk) #18679   2

#[Counts] cor between sc and bulk; Not good
cor.test(countBulk$bulkP1, countBulk$bulkP2) # 0.9943626
plot(countBulk$bulkP1, countBulk$bulkP2)


# [cpm] cor between sc and bulk
getNormalizedCts <- function ( cts ) {
  #cts <- read.table ( ctsPath , header = T , as.is = T )
  as.data.frame( apply ( cts , 2 , function ( x ) { log2 ( ( 10^6 ) * x / sum ( x ) + 1 ) }) )
}

cpmBulk=getNormalizedCts(bulkAll0)
head(cpmBulk)
dim(cpmBulk) #57792     2
#
cpmBulk2=cpmBulk[genes,]
dim(cpmBulk2) #18679     2

##
dim(sc0)
sc.cpm=as.data.frame( apply(sc0, 2, function(x){
  x/sum(x)*1e6
}) )
sc.cpm.Pool=as.data.frame( apply(sc.cpm, 1, sum) )
colnames(sc.cpm.Pool)=c('Pool')
sc.cpm.Pool$gene=rownames(sc.cpm.Pool)
dim(sc.cpm.Pool) #18679     2
head(sc.cpm.Pool)
#
cpmScPool2=data.frame(
  gene=sc.cpm.Pool$gene,
  Pooled=log2(1+sc.cpm.Pool$Pool/255)
)
row.names(cpmScPool2)=cpmScPool2$gene
head(cpmScPool2)

## cor of bulk1 and bulk2
cor.test(cpmBulk2$bulkP1, cpmBulk2$bulkP2) # 0.9819343 

## cor of bulk1 or bulk2 to scPooled
cor.test(cpmBulk2$bulkP1, cpmScPool2$Pooled) #0.9573284 
cor.test(cpmBulk2$bulkP2, cpmScPool2$Pooled) #0.9576741

plot(cpmBulk2$bulkP1, cpmBulk2$bulkP2)
plot(cpmBulk2$bulkP1, cpmScPool2$Pooled)


#########
cor0=cor(cpmBulk2$bulkP1, cpmBulk2$bulkP2); print(cor0) #0.981
cor1=cor(cpmBulk2$bulkP1, cpmScPool2$Pooled); print(cor1) #0.957

library(ggplot2)
g=ggplot(NULL, aes(cpmBulk2$bulkP1, cpmScPool2$Pooled ) )+theme_classic()+
  geom_point(size=0.02, alpha=0.1)+
  geom_text(aes(x=0,y=14,label=paste0("r=", round(cor1[1],3) ) ),hjust=0 )+  
  labs(x="Tag of bulk data (log2)", y="Tag of pooled single cells (log2)",
       title="Gene expresion in log2(cpm+1)")
g

g2=ggplot(NULL, aes(cpmBulk2$bulkP1, cpmBulk2$bulkP2 ) )+theme_classic()+
  geom_point(size=0.02, alpha=0.1)+
  geom_text(aes(x=0,y=14,label=paste0("Rp=", round(cor0[1],3) ) ),hjust=0 )+  
  labs(x="MDA-MB-468 bulk1", y="MDA-MB-468 bulk2",
       title="Gene expresion in log2(cpm+1)")
g2


CairoPDF('bulk_PoolSC.cor-2-new.pdf', width=3.3,height=3.3)
print(g)
print(g2)
dev.off()

