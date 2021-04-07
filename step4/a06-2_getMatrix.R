#aim: combine small bed to get APA matrix( very slow! R version is too slow! use shell!)
# v0.1 can run
# v0.2 correct: id=id_list[1]
# v0.3 for cell line

setwd('/home/wangjl/data/apa/190705PAS/bed/freq2/')
getwd()


# read cell id list
id_list=readLines('/home/wangjl/data/apa/190705PAS/225.cellID')
#id_list=id_list[1:10] #debug

head(id_list) # "c01_ROW07" "c01_ROW12"
length(id_list) #225


# get first bed file, as a frame, other file load 2nd column
main_matrix=read.table(paste0(id_list[1],".freq"),header=F)
main_matrix[1:4,]
cb=id_list[1]
colnames(main_matrix)=c('PAS',cb)
head(main_matrix)


#function: get 2nd col of each cell
getInter=function(id){
  info=read.table(paste0(id,".freq"),header=F)
  info[,2]
}

# add all other col
for( i in 2:length(id_list) ){
  if(i %% 2==0) print(i)
  #if(i>100) break;
  cb=id_list[i];
  cbdata=getInter(cb);
  main_matrix[,cb]=cbdata
}

dim(main_matrix) #[1] 1559047      11

# write to file, sep=tab
write.table(main_matrix,"all225_matrix_APACountsV2.txt",row.names=F,sep="\t")

#end
print("====End===")




##QC0:compare result with shell
sum(main_matrix)
dim(main_matrix)
main_matrix[1:3,1:6]

wjlR=read.table('/data/jinwf/wangjl/apa/190705PAS/bed/freq2/all225_matrix_APACountsV2.txt', header = T,row.names = 1)
sum(wjlR)
dim(wjlR)
wjlR[1:3,1:6]

wjl=read.table('/data/jinwf/wangjl/apa/190705PAS/bed/freq2/test2.txt', header = T,row.names = 1)
sum(wjl)
dim(wjl)
wjl[1:3,1:6]

#
wjl2=wjlR-wjl
sum(wjl2)
# the same
rm(wjl,wjl2,wjlR)


## QC1: at least 2 cells with at least 1 reads support
apa=read.table('/data/jinwf/wangjl/apa/190705PAS/bed/freq2/all225_matrix_APACountsV2.txt', header = T,row.names = 1)
dim(apa) #[1] 742288    225
#
table(apply(apa>1,1,sum)>1) #[1] 742288F    225T
table(apply(apa>0,1,sum)>1) #[1] 386653 355635
#
apa2=apa[apply(apa>1,1,sum)>1,]
dim(apa2) #[1] 267961    225
sum(apa) #[1] 224548724
sum(apa2) #[1] 215435060
sum(apa2)/sum(apa)*100 #95.94%

# write to file, sep=tab
write.table(apa2,"all225_matrix_APACountsV3.txt",sep="\t")

# save APA site col
pas=rownames(apa2)
head(pas)
length(pas) #267961
writeLines(pas,'all225_pas_267961.txt')


## QC2: no genomic 4A after polyA site
dim(apa2) #267961    225
apa[1:4,1:4]
# read PAS list
keep=read.table('/home/wangjl/data/apa/190705PAS/bed/pasPosition_no4A.bed',header = F,sep="\t")
dim(keep) #[1] 440760      6
keep[1:4,]
keep=keep$V4
keep=intersect(keep, rownames(apa2))
length(keep)
apa4=apa2[keep,]
dim(apa4) #[1] 150626    225

# write to file, sep=tab
write.table(apa2,"all225_matrix_APACountsV4.txt",sep="\t")
pas=rownames(apa4)
head(pas)
length(pas) #150626
writeLines(pas,'all225_pasV4_150626.txt')


# check
sum(apa4['chr12:6647519:+',])


# plot counts per Site, in log scale
readPerPas=apply(apa4,1,sum)
length(readPerPas) #[1] 150626
readPerPas[1:20]
#
ggplot(NULL, aes(x=log10(unname(readPerPas) ) ))+#geom_density()+
  geom_histogram(binwidth = 0.03,fill="white",colour="red")+
  xlim(0,5.5)+
  geom_vline(aes(xintercept=1.2),linetype="dashed")+
  labs(x='log(Reads per polyA site)')
#
qplot(log10(unname(readPerPas)), geom="histogram", binwidth=0.5)

#
dt=unname(readPerPas)
table(unname(readPerPas)>0)   #        150626T
table(unname(readPerPas)>10)  #F19604  131022T
table(unname(readPerPas)>100) #F78994   71632T
table(unname(readPerPas)>1e3) #F136723  13903T
table(unname(readPerPas)>1e4) #F148420   2206T
table(unname(readPerPas)>1e5) #F150430    196T
table(unname(readPerPas)>1e6) #F150618      8T
table(unname(readPerPas)>1e7) #F150626
readPerPas[dt>1e6]

#10**1.2=15.8489
table(unname(readPerPas)>10**1.2) #F28480 122146T


# need further filter? Risk: rising filter may filter out low expressed TF
table(apply(apa4>1,1,sum)>1) #       150626T 
table(apply(apa4>1,1,sum)>2) #F 58256 92370T
table(apply(apa4>1,1,sum)>3) #F 82688 67938T
table(apply(apa4>1,1,sum)>4) #F 95931 54695T # cell>=5, read per cell >=2
table(apply(apa4>1,1,sum)>5) #F104520 46106T
table(apply(apa4>1,1,sum)>10)#F122958  27668T
#225*10%=22.5
table(apply(apa4>1,1,sum)>22) #F134986  15640T

#
# read per cell, not influence so much
table(apply(apa4>2,1,sum)>4) #F105093  45533T
table(apply(apa4>3,1,sum)>4) #F109867  40759T

#
ggplot(NULL, aes(x=log10(unname(apply(apa4,1,sum)) ) ))+#geom_density()+
  geom_histogram(binwidth = 0.03,fill="white",colour="red")+
  xlim(0,5.5)+
  geom_vline(aes(xintercept=1.2),linetype="dashed")+
  labs(x='log10(Reads per polyA site)')
#



#QC1.5 enhence: filter again 
keep2=apply(apa4>=3,1,sum)>=5;table(keep2) # (cell>=5, read per cell >=3)
#105093  45533
apa5=apa4[keep2,]
ggplot(NULL, aes(x=log10(unname(apply(apa5,1,sum)) ) ))+#geom_density()+
  geom_histogram(binwidth = 0.03,fill="white",colour="red")+
  xlim(1,6)+
  #geom_vline(aes(xintercept=1.2),linetype="dashed")+
  labs(x='log10(Reads per polyA site)', title='PAS supported by at 3 reads per cell in at least 5 cells(45533PAS)')
#

# Rising filter shreshold, may losing important TF, but the result maybe more reliable.
pas=rownames(apa5)
head(pas)
length(pas) #45533
writeLines(pas,'all225_pasV5_45533.txt')#filter too much?!

#sum(apa['chr10:76991020:+',])
