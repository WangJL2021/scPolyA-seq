# get no correlated genes
setwd("/data/jinwf/wangjl/apa/20200701Fig/f6/result/cor_BC/")
getwd()

# load data
cor_df=read.table("MBA-MD-468_01-volcano_Spearman_cor_gDPAU_RNA.df.txt")
cor_df=cor_df[order(-cor_df$cor),]
head(cor_df)
table(cor_df$sig)
#n.s. negative positive 
#3181      595      222 

table(cor_df$sig, cor_df$cor>0)
#          FALSE TRUE
#n.s.      1945 1236
#negative   595    0
#positive     0  222






############
# all n.s. genes
gene.ns=rownames(cor_df[which(cor_df$sig=="n.s."),])
writeLines( gene.ns, "MBA-MD-468_02-Spearman_cor_gDPAU_RNA.No-Correlation.gene.txt")

# plus
writeLines( rownames(cor_df[which(cor_df$sig=="n.s." & cor_df$cor>0),]), "no-cor-plus.gene.txt")

#minus
writeLines( rownames(cor_df[which(cor_df$sig=="n.s." & cor_df$cor<0),]), "no-cor-minus.gene.txt")


# hist overlay
hist(cor_df$cor, n=150, col=rgb(0,0,1,0.2), main="Histgram of correlation")
hist(cor_df[gene.ns, ]$cor, n=150, col=rgb(1,0,0,0.2), add=T)

# try again
hist( log10( 1e-30+ abs(cor_df$cor)), n=100, col=rgb(0,0,1,0.2), 
      main="Histgram of correlation")
hist( log10( 1e-30+ abs(cor_df[gene.ns, ]$cor)), n=100, col=rgb(1,0,0,0.2), add=T)



############
# all n.s. genes in to bin
hist(cor_df[gene.ns, ]$cor, n=150, col=rgb(1,0,0,0.2))
# top 5% genes
df2=cor_df[gene.ns, ]
df2$cor2=abs(df2$cor)
df2=df2[order(-df2$cor2),]
head(df2)
dim(df2) #3181 7
#
hist(df2$cor, n=150, col=rgb(1,0,0,0.2)) #-0.4 to 0.4
hist(df2$cor2, n=150, col=rgb(1,0,0,0.2))
dim(df2) #3181 7

n=round( nrow(df2[gene.ns, ]) * 0.1); n #159
#
writeLines( rownames(df2)[1:n], "noSig_cor-abs_High.gene.txt")
writeLines( rownames(df2)[(3181-n):3181], "noSig_cor-abs_Low.gene.txt")
#
# by bin
writeLines( rownames(df2[which(df2$cor2<0.05),]), "noSig_cor-abs_lt0.05.gene.txt")

writeLines( rownames(df2[which(df2$cor2<0.2),]), "noSig_cor-abs_lt0.2.gene.txt")
writeLines( rownames(df2[which(df2$cor2>=0.2),]), "noSig_cor-abs_gt0.2.gene.txt")
#
writeLines( rownames(df2[which(df2$cor2>=0.3),]), "noSig_cor-abs_gt0.3.gene.txt")

#
table(df2$cor2<0.02)
writeLines( rownames(df2[which(df2$cor2<0.01),]), "noSig_cor-abs_lt0.01.gene.txt") #153
writeLines( rownames(df2[which(df2$cor2<0.02),]), "noSig_cor-abs_lt0.02.gene.txt") #279

##############
# sig neg
df3=cor_df[which(cor_df$sig=="negative"),]
dim(df3) #595
hist(df3$cor, n=100)
#
table( df3$cor>=(-0.45) )
#FALSE  TRUE 
# 208   387 
writeLines( rownames(df3[which(df3$cor<(-0.45) ),]), "negSig_cor-lt-0.45.gene.txt")
writeLines( rownames(df3[which(df3$cor>=(-0.45) ),]), "negSig_cor-gt-0.45.gene.txt")


         
