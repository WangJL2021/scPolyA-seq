# aim: get expression from 293T cells data
#

setwd("/home/wangjl/data/apa/20200701Fig/f1/3T3exp/")
getwd()
#
library(dplyr)
library(Seurat)
library(patchwork)

# Load the scRNA dataset
scRNA.data <- Read10X(data.dir = "/home/wangjl/data/ref/293T/filtered_matrices_mex/hg19/") #
#scRNA.data <- Read10X(data.dir ="/data/jinwf/wangjl/ref/pbmc10k/filtered_feature_bc_matrix") #v3.1 cellrange4.0.0
#
dim(scRNA.data) #32738  2885
head(scRNA.data)
n=1;scRNA.data[1:100, n:(40+n)]


# Initialize the Seurat object with the raw (non-normalized data).
scRNA <- CreateSeuratObject(counts = scRNA.data, project = "293T", min.cells = 3, min.features = 200)
dim(scRNA@assays$RNA@counts) #16316  2885

hist(apply(scRNA@assays$RNA@counts, 2, sum), n=100) #counts
hist(apply(scRNA@assays$RNA@counts>0, 2, sum), n=100) #gene number

# write to file
write.table(scRNA@assays$RNA@counts, '293T_exp.dt.txt')




########
#(1) plot counts / gene per gene
rnaM=scRNA@assays$RNA@counts
dim(rnaM) #16316  2885
countsPerCell=apply(rnaM, 2, sum)
genePerCell=apply(rnaM>0, 2, sum)


library(ggplot2)
ggplot(data.frame(value=countsPerCell), aes(x=1, y=value/1e6))+geom_boxplot()+
  geom_jitter(alpha=0.3)+
  labs(x="293T cells(by 10x)", y='Million counts per cell')+
  theme_classic()
#

ggplot(data.frame(value=genePerCell), aes(x=1, y=value))+geom_boxplot()+
  geom_jitter(alpha=0.3)+
  labs(x="293T cells(by 10x)", y='Gene number per cell')+
  theme_classic()
#

## 
# read RNA matrix of C1
rnaMc1=read.csv('/data/jinwf/wangjl/apa/20200701Fig/f2/BC_HeLa/BC_HeLa.222cells.count.V4.csv', row.names = 1)
head(rnaMc1[,1:10])
dim(rnaMc1) #18662   222

#
countsPerCell2=apply(rnaMc1, 2, sum)
genePerCell2=apply(rnaMc1>0, 2, sum)
#
hist(countsPerCell2, n=100)
hist(genePerCell2, n=100)


####################
# counts
df1=data.frame(
  tech=c( rep('10x', length(countsPerCell)) ,   rep('c1', length(countsPerCell2))  ),
  value=c(as.numeric(countsPerCell),  as.numeric(countsPerCell2))
)
g1=ggplot(df1, aes(x=tech, y=value/1e6, color=tech))+geom_boxplot(outlier.alpha = 0)+
  #geom_jitter(size=0.1, alpha=0.2)+
  scale_x_discrete(labels = c('10x scRNA-seq', 'scPolyA-seq')) +
  labs(x="", y='Million counts per cell')+
  theme_classic()+theme(legend.position="none")+
  theme(axis.text.x=element_text(angle=60, hjust=1,size=10) ); g1
#
#
# gene
df2=data.frame(
  tech=c( rep('10x', length(genePerCell)) ,   rep('c1', length(genePerCell2))  ),
  value=c(as.numeric(genePerCell),  as.numeric(genePerCell2))
)
g2=ggplot(df2, aes(x=tech, y=value, color=tech))+geom_boxplot(outlier.alpha = 0)+
  #geom_jitter(size=0.1, alpha=0.05)+
  scale_x_discrete(labels = c('10x scRNA-seq', 'scPolyA-seq')) +
  labs(x="", y='Gene number per cell')+
  theme_classic()+theme(legend.position="none")+
  theme(axis.text.x=element_text(angle=60, hjust=1,size=10) ); g2
#

# combine to one
library(gridExtra)

CairoPDF('10x_c1_gene_counts_per_cell-2.pdf', width=3, height=3)
grid.arrange(g1,g2, nrow=1)
dev.off()
###

