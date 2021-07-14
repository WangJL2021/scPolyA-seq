input=snakemake@input[[1]]
output=snakemake@output[[1]]

# get data frame
df=read.table( input, header=T)
df$ratio=df$uniqReads/df$total
df$multiple=df$total - df$uniqReads
# order by total reads number
df=df[order(df$total),]
head( df )

# fig 1: barpot
colNames=c("uniqReads","multiple");
mat= as.matrix( df[, colNames] )
cols=c('purple', 'grey')

pdf(output, width=4, height=3.5)
oPar=par(); par(mar=c(4,4,2,6))
barplot( t(mat), col=cols, border=NA, ylab="Reads per cell" )
xy=par('usr')
legend(x=xy[2], y=xy[4], fill=rev(cols), legend= rev(colNames), 
       border = NA,
       bty="n",
       xpd=T )

# fig 2: plot
par(oPar)
plot(df$ratio, type='o', lwd=2,
     xlab="Cell id", ylab="Uniq map ratio", 
     #ylim=c(y1,y2),
     main="Uniq mapped ratio per cell")
dev.off()

