#visual: counts per cell after split by columns.
#draw reads per well, in log scale
#v0.1 heatmap in log for color, label million read in each well.
#v0.2 output in pdf


#############################
#input data
setwd("/home/wangjl/data/apa/fq_files")
d1=read.csv("demultiplex_report.xls",sep="\t")
head(d1);dim(d1) #492   3
d2=read.csv("demultiplex_report2.xls",sep="\t")
head(d2);dim(d2) #164   3
#
d3=rbind(d1,d2)
head(d3);dim(d3) #[1] 656   3
rm(d1,d2)

aa=d3
######
hist(aa$read_number,n=492,xlab="Read Counts per Cell",main="")
# outliers
aa[which(aa$read_number>1e7),] 
#     col_sample_name cell_sample_name read_number
#126             c15        c15_ROW02    18328816
#444              c8         c8_ROW33    10257039
#645             c10        c10_ROW29    12108205
#648             c10        c10_ROW32    14913383
#653             c10        c10_ROW37    10568993


library("Cairo")
CairoPDF(file='/data/jinwf/wangjl/apa/190515/whatIs_EmptyWell.pdf',width=6,height=6)
par(mar=c(4.5, 4, 0.5, 1.5))
par(mfrow=c(2,1))
#fig1 
hist(aa$read_number,n=492,xlab="Read Counts per Cell",ylab="Freq",main="",
     xlim=c(0,2e7))
#draw lines
abline(v=10**5.5,col='blue',lty=2)
abline(v=5.5e6,col='red',lty=2)
#text
text(-4e5,20,"white",cex=0.8)
text(2.5e6,20,"single",cex=1,col="blue")
text(1e7,20,"dual?",cex=1.2,col="#D42331")
#
text(5.4e6,25,"x=5.5e6",cex=0.7,col="red",pos=4)
#end

#fig2 draw again in log scale 
hist(log10(aa$read_number),n=500, xlab="log10[Read Counts per Cell]",ylab="Freq",main="")
#
abline(v=5.5,col='blue',lty=2)
abline(v=log10(5.5e6),col='red',lty=2)
#
text(5,7,"white",cex=0.8)
text(6,7,"single",cex=1,col="blue")
text(7.1,7,"dual?",cex=1.2,col="#D42331")
#
text(5.47,4,"x=1e5.5",cex=0.7,col="blue",pos=4)
#end2
dev.off()



# get row number
aa2=aa;
aa2$cell_sample_name=gsub('-','_',aa2$cell_sample_name) #-替换为_
d1=c()
for(i in 1:nrow(aa2)){
  x=as.character(aa2[i,2])
  #d1=c(d1,substr(x,nchar(x)-1,nchar(x)) )
  d1=c(d1, strsplit(x,"_")[[1]][2] )
}
aa2$row=d1
aa3=aa2[,c(1,4,3)]
head(aa2)
head(aa3)


#long2wide
library(reshape2)
aa3=dcast(aa3, row~col_sample_name,value.var="read_number")
head(aa3)
# the 1st columns is rowname, del it
aa4=aa3[,2:ncol(aa3)]
head(aa4)
rownames(aa4)=aa3$row
head(aa4)
class(aa4)
#re order
aa5=aa4[rev(rownames(aa4)),c("c01","c2","c3","c04","c05","c6","c7","c8","c9","c10",
                             "c12","c13","c14","c15","c16","c19")]
head(aa5)
# rm(aa2,aa3,aa4)
dim(aa5) #41 16


#############################
#visual
data=t(aa5)
 
m1=log10(min(data))
m2=log10(max(data))
breaks.frequency = seq(from=m1, to=m2, length.out=20)
myColors = colorRampPalette(c("white", "#6CBF4B","#D42331"))

CairoPDF(file="/data/jinwf/wangjl/apa/190515/heatmap_eachWell_by_logCounts.pdf",width=8,height = 8);
par(mar=c(5,5,4,2))
image(1:nrow(data), 1:ncol(data), as.matrix(log10(data)),  
      breaks=breaks.frequency,#color+1, increase
      col=myColors(length(breaks.frequency)-1), 
      axes = F,
      cex = 0.5,
      main="log10(Million Reads per Well)",
      xlab = "", ylab = "")
#
axis(3, #top of x axis
     mgp=c(1,0.2,0), 
     at=1:nrow(data), 
     rownames(data), 
     las=2,
     cex.axis=0.6,
     lwd=0,lwd.ticks=0) #
axis(2, at=1:ncol(data), colnames(data), 
     mgp=c(1,0.2,0),
     las=2, cex.axis=0.6,
     lwd=0,lwd.ticks=0)
# heatmap: text on the well
for (x in 1:nrow(data)) {
  for (y in 1:ncol(data))  {
    text(x, y, round(data[x, y]/1e6,2),cex=0.6) 
  }
}
dev.off()
#

## from the image, the threshold of empty is 0.23e6
head(aa);dim(aa)


################################
#cells filtered out: cellID
df.empty=aa[aa$read_number<0.23e6,]
dim(df.empty);head(df.empty) #[1] 70   3

#Left
df.emptyLeft=df.empty[df.empty$col_sample_name %in% c("c01","c2","c3","c04","c05",
                                         "c6","c7","c8","c9","c10"),]
dim(df.emptyLeft);head(df.emptyLeft) #37 3
df.emptyLeft
df.emptyLeft$cell_sample_name
#
#Right
df.emptyRight=df.empty[df.empty$col_sample_name %in% c("c12","c13","c14","c15",
                                                       "c16","c19"),]
dim(df.emptyRight);head(df.emptyRight) #33 3
df.emptyRight
df.emptyRight$cell_sample_name
#
#save
write.table(df.emptyLeft$cell_sample_name,file="../190515/01empty.Left.cellID",
          quote=F,row.names = F,col.names = F)
write.table(df.emptyRight$cell_sample_name,file="../190515/01empty.Right.cellID",
            quote=F,row.names = F,col.names = F)
#



################################
#cellID that keep
df=aa[aa$read_number>=0.23e6,]
dim(df);head(df) #[1] 586   3
#去掉Undetermined行
df=df[grep("Undetermined",df$cell_sample_name,invert=T),]
dim(df)#570 3

#Left
df.Left=df[df$col_sample_name %in% c("c01","c2","c3","c04","c05",
                                                      "c6","c7","c8","c9","c10"),]
dim(df.Left) #363
#Right
df.Right=df[df$col_sample_name %in% c("c12","c13","c14","c15",
                                                       "c16","c19"),]
dim(df.Right) #207
#save
write.table(df.Left$cell_sample_name,file="../190515/01NO_empty.Left.cellID",
            quote=F,row.names = F,col.names = F)
write.table(df.Right$cell_sample_name,file="../190515/01NO_empty.Right.cellID",
            quote=F,row.names = F,col.names = F)
