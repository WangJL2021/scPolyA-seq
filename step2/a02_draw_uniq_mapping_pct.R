# map to hg19 / mm10, uniq counts plot

setwd("/home/wangjl/data/apa/190515")
getwd()
#source("https://bioc.ism.ac.jp/biocLite.R")
#biocLite("ggplot2")
library("ggplot2")
library(Cairo)

###########################
#load data
d1=read.table("hg19_mm10_uniqCount.log", sep="\t")
head(d1);dim(d1)

df=NULL
for(i in 1:nrow(d1)){
  x=as.character(d1[i,2]);x
  #get number before keyword
  getNumBefore=function(kw,xx=x){
    x1=strsplit(xx,kw)[[1]][1];x1
    x2=strsplit(x1," ")[[1]];x2
    as.numeric(x2[length(x2)]);
  }
  hg19=getNumBefore("HUMAN")
  #
  mm10=getNumBefore("MOUSE")
  df=rbind(df, data.frame(
    cid=as.character(d1[i,1]),
    hg19=hg19,
    mm10=mm10
  ))
}
head(df)
dim(df) #400 3
#
# get pct data
df$hg19p=df$hg19/(df$hg19+df$mm10)
df$mm10p=df$mm10/(df$hg19+df$mm10)
head(df)


###########################
# draw scatter plot
#boxplot(df$hg19,df$mm10) #not ok
#hist(df$hg19/df$mm10, n=400) #not ok
write.table(df, "pic00_uniqCounts_hg19_mm10-plot.df.txt")

CairoPDF(file="pic00_uniqCounts_hg19_mm10-plot.pdf",width=3,height=3)
# normal scale
ggplot(df, aes(mm10/1e6,hg19/1e6))+geom_point(alpha=0.2)+
  theme_bw()+
  labs(title="Unique mapping counts[rm empty]")
  #,x="Counts on mouse genome",y="Counts on human genome"

#log scale
ggplot(df, aes(log10(mm10/1e6),log10(hg19/1e6) ))+geom_point(alpha=0.2)+
  theme_bw()+
  labs(title="Unique mapping counts(log scale)[rm empty]",
       x="log10[mm10 counts]",y="log10[hg19 Counts]")

#pct relative(not much info)
ggplot(df, aes(mm10/(mm10+hg19)*100,hg19/(mm10+hg19)*100))+
  geom_point(color="red",alpha=0.2)+
  geom_rug()+
  xlim(0,100)+ylim(0,100)+
  theme_bw()+
  labs(title="Unique mapped counts ratio[rm empty]",x="mm10",y="hg19")
dev.off()






# test 0.8, 0.9
df__=df
df__$source="mixed"
df__[which(df__$hg19p>0.9),]$source="hg19"
df__[which(df__$mm10p>0.9),]$source="mm10"
table(df__$source)
# hg19 mixed  mm10 
# 13    35   352

df$source="mixed"
df[which(df$hg19p>0.8),]$source="hg19"
df[which(df$mm10p>0.8),]$source="mm10"
table(df$source)
#hg19 mixed  mm10 
#13    20   367 
table(df__$source, df$source)
#      hg19 mixed mm10
#hg19    13     0    0
#mixed    0    20   15
#mm10     0     0  352
# use 0.8 later.

# remove empty well, on the left side.
cid=readLines("01NO_empty.Left.cellID")
df2=df[which(df$cid %in% cid),]
tb=table(df2$source);tb
#hg19 mixed  mm10 
#13    15   335 
g=ggplot(df2, aes(mm10/1e6,hg19/1e6))+geom_point(alpha=0.4)+
  labs(title="Unique mapping counts")+
  labs(x="Million counts on mm10", y="Million counts on hg19")+
  theme_bw();g

g1=ggplot(df2, aes(mm10/1e6,hg19/1e6,color=source))+geom_point(alpha=0.4)+
  labs(title="Unique mapping counts")+
  labs(x="Million counts on mm10", y="Million counts on hg19")+
  theme_bw()+
  scale_color_discrete(name  ="",breaks=c("hg19", "mixed", "mm10"),
                       labels=c("hg19(13)", "mixed(15)", "mm10(335)"));g1

## save
write.csv(df,  'hg19_mm10_plot/uniqCounts_hg19_mm10-plot-df.csv')
write.csv(df2, 'hg19_mm10_plot/uniqCounts_hg19_mm10-plot-df2.csv')

CairoPDF(file="hg19_mm10_plot/uniqCounts_hg19_mm10-plot-v1.pdf",width=2.4,height=2.5)
print(g)
dev.off()

CairoPDF(file="hg19_mm10_plot/uniqCounts_hg19_mm10-plot-v2.pdf",width=3.5,height=2.5)
print(g1)
print(g1+xlim(0,12)+ylim(0, 12))
dev.off()










# melt to one column
library(reshape2)
d1=melt(df[,c(1,2,3)],id.vars = "cid", neasure.vars=c("hg19",'mm10'),
        variable.name="Genome",value.name="counts")
head(d1)

#get pct data
d2=melt(df[,c(1,4,5)],id.vars="cid", neasure.vars=c("hg19p",'mm10p'), 
        variable.name="Genome",value.name="pct")
head(d2)

#test
#ggplot(d1, aes(x=cid, y=counts,fill=Genome))+geom_bar(stat="identity")

#uniq mapped reads, total counts
CairoPDF(file="pic01_uniqCounts_distribution.pdf",width=6,height=3)
#tmp2=tmp[order( tmp$log.P., tmp$Term),]
ggplot(d1, aes(x=reorder(cid,counts), y=counts/1e6,fill=Genome))+
  geom_bar(stat="identity")+
  labs(title="",x="Cell ID",y="Uniq mapped counts(Million)")+
  scale_x_discrete( labels=NULL)
#
ggplot(d1, aes(x=reorder(cid,counts, function(x){x[2]}), y=counts/1e6,fill=Genome))+
  geom_bar(stat="identity")+
  labs(title="",x="Cell ID",y="Uniq mapped counts(Million)")+
  scale_x_discrete(labels=NULL)
dev.off()


#order 
tmp=df[order( -df$hg19p, -df$hg19),]
head(tmp)

#uniq mapped reads, order by mm10
CairoPDF(file="pic02_uniqCounts_distribution-bymm10.pdf",width=6,height=3)
ggplot(d1, aes(x=reorder(cid,counts, function(x){x[2]}), y=counts/1e6,fill=Genome))+
  geom_bar(stat="identity")+
  labs(title="",x="Cell ID",y="Uniq mapped counts(Million)")+
  scale_x_discrete(limits=tmp$cid, labels=NULL)
#dev.off()
#uniq mapped reads: pct, oder by hg19
#CairoPDF(file="pic01_uniqCounts_distribution-bymm10_100pct.pdf",width=7,height=3)
ggplot(d2, aes(x=reorder(cid,pct, function(x){x[2]}), y=pct*100,fill=Genome))+
  geom_bar(stat="identity", position="fill")+
  labs(title="",x="Cell ID",y="Uniq mapped counts ratio")+
  scale_x_discrete(limits=tmp$cid, labels=NULL)
dev.off()


# test for hg
df2=df[df$hg19p>0.9,];dim(df2) #[1] 13  5
df2

#save
save(df,tmp, file='uniqMappedCounts_before_rmdup.RDS')




############
load("uniqMappedCounts_before_rmdup.RDS")
head(df);dim(df)
df$total=df$hg19+df$mm10
df2=df[order(df$total),]
head(df2)
df3=df2[df2$total>=0.22e6,]
head(df3);dim(df3) #359 6
#
df.empty=df[df$total<0.22e6,]
head(df.empty);dim(df.empty) #41 6
df.empty
