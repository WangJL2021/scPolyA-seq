#b01 filter out cells with too many counts

setwd("/data/jinwf/wangjl/apa/190517R/")
d1=read.table("uniq_mapped_reads.207.txt",header=F)
head(d1);dim(d1) #207 6

df=data.frame(
  rname=d1$V1,
  counts=d1$V3,
  pct=d1$V5,
  row.names = "rname"
)
head(df);dim(df) #207 2
df2=df[order(df$counts),]

#plot
library("Cairo")
CairoPDF(file="b01_Uniquely mapped reads number.pdf", width=10,height=10)
par(mar=c(5,4,4,8)+0.1)
tmp2=barplot(df2$counts/1e6,xlab="Cell ID",ylab="Million Reads",
        main="Uniquely mapped reads number[bc+HeLa]")
#threshold
th1=median(df2$counts)-mad(df2$counts)*3;th1
th2=median(df2$counts)+mad(df2$counts)*3;th2
#
abline(h=th1/1e6,lty=2,col="red")
abline(h=th2/1e6,lty=2,col="red")
#
df3=df2
df3$pct=as.numeric( gsub("%","",df3$pct) )
df3$pctA=(df3$pct-min(df3$pct) )/( max(df3$pct)-min(df3$pct) ) #76.94-65.79
#
points(tmp2,df3$pctA*4+6, type="l",col="blue",lwd=0.1)
kd=seq(min(df3$pct), max(df3$pct), ( max(df3$pct)-min(df3$pct) )/4)
axis(4,at=seq(6,10,1),labels=round(kd,2),
     col="blue",col.axis="blue",las=2,cex=0.7,tck=-0.02)
dev.off()



#cell id
tmp=df2[df2$counts>th2,];tmp
#            counts    pct
#c12_ROW32  5383983 75.08%
#c14_ROW04  5971786 72.42%
#c14_ROW03  6668614 73.63%
#c15_ROW02 12276181 69.39%

#rm posible doublet
write.table(tmp,file="02FilteredByUniqMapCounts-discard.4.Right.csv",
            quote=F,row.names = T,col.names = T)
write.table(rownames(df2[df2$counts<=th2,]),file="02FilteredByUniqMapCounts.203.Right.cellID",
            quote=F,row.names = F,col.names = F)
