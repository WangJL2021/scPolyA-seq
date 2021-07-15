#
setwd("/data/jinwf/wangjl/apa/20200701Fig/f4/cycles/")
getwd()

output="sync_normal/"
keyword="sync_vs_normal"
# load data ====
df1_all=read.table("sync_normal/df1_gDPAU_delta.df.txt")
head(df1_all)

df2_all=read.table("sync_normal/df2_gDPAU_chisqP.df.txt")
head(df2_all)


# combine df ====
gene.common=intersect(rownames(df1_all), rownames(df2_all))
df1=df1_all[gene.common, ]
df2=df2_all[gene.common, ]
dim(df1) #3940 12
dim(df2)
colnames(df2)=c("gene",    "countsA", "countsB", "p.Chisq",       "padj.Chisq")

df=cbind(df1, df2)
head(df)

# plot =====
df$sig="ns"
df[which( (df$delta)>30  & df$padj.Chisq<0.05), ]$sig="up"
df[which( (df$delta)<(-30)  & df$padj.Chisq<0.05), ]$sig="down"
table(df$sig)
#down   ns   up 
#203 3526  211 
# save
write.table(df, "sync_normal/df_gDPAU_sync_normal_chisqP_delta.df.txt")


tb=table(df$sig); print(tb)

CairoPDF(paste0(output, '02_volcano_gDPAU_chisqP_',keyword,'-2.pdf'),width=3.3,height=3.5)
g=ggplot(df, aes(delta/100, -log10(padj.Chisq), color=factor(sig, levels=c('down', 'ns', 'up')) ) )+
  geom_point(size=0.1, alpha=0.5)+theme_bw()+
  labs( title=paste("gDPAU:", keyword, length(cid.A),length(cid.B)), 
        x='Change of gDPAU', y="-log10(adj.p)" ) + theme_bw()+
  theme(legend.box = "horizontal",
        legend.key.size=unit(6,"pt"),
        legend.position="bottom") +
  scale_color_manual('', labels=c( paste0("shorten(",tb[1],")"),
                                   paste0('n.s.(',tb[2],')'),paste0("lengthen(",tb[3],")") ), 
                     values=c("blue", "#dddddd",'red') );
dd_text=rbind(
  (function(){
    df2=df[which(df$sig=='up'),]
    head(df2[order(df2$padj), ])
  })(),
  (function(){
    df2=df[which(df$sig=='down'),]
    head(df2[order(df2$padj), ])
  })()
);pd(dd_text)
g2=g+geom_text_repel(data=dd_text, aes(x=delta/100, y=-log10(padj.Chisq), label=rownames(dd_text)),
                     color="black",size=3,alpha=0.6)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=2)))
print(g2)
dev.off()
