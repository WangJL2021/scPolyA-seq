# Code for scPolyA-seq

> scPolyA-seq produced a strand-specific, 3’end RNA-seq library.




## scPolyA-pipe

The pipeline to get polyA site and single-cell polyA site matrix from raw fastq of scPolyA-seq.


**How to use scPolyA-pipe to do *de-novo* polyA site discovery and get APA matrix from scPolyA-seq?**


### 1. Firt download this set of snakemake rules to anywhere on your Linux system, 
and install some required softwares.

```
Our tests based on the following envirenments:
OS: CentOS 7
$ python3 -V #Python 3.7.3
$ snakemake --version  #6.4.1
$ R --version  #R version 3.6.0 (2019-04-26) -- "Planting of a Tree"
$ STAR --version #STAR_2.5.2b
$ fastqc --version  #FastQC v0.11.8
$ multiqc --version #multiqc, version 1.7

$ samtools --version
samtools 1.9-69-gb217a91
Using htslib 1.9-149-gf5b75ff

$ bedtools bamtobed --version
Tool:    bedtools bamtobed (aka bamToBed)
Version: v2.25.0
```



### 2. Goto your working directory, put your fastq file

```
for C1, put your cellId_R2.fastq.gz in raw/
raw/c01ROW07.fastq.gz
raw/c01ROW12.fastq.gz

or fastq file in fq/
fq/c01ROW07.fq
fq/c01ROW12.fq
```



### 3. Make a new file named config.yaml under working directory, 
put important initial settings in this config file.



```
$ head config.yaml
Snk_RootDir: /home/wangjl/data/scPolyA-seq/scPolyA-pipe
star_index: /home/wangjl/data/ref/hg19/index/star/

RawPath: /home/wangjl/data/scPolyA-seq/data/raw_sc/human
test: /home/wangjl/data/apa/test/raw
Samples:
 - c01ROW07
 - c01ROW12
```


### 4. (We recommend open an Tmux panal, in case the pipeline ceased on network error) 
Run the pipeline in the working directory:

```
$ snakemake -s /home/wangjl/data/scPolyA-seq/scPolyA-pipe/main.sf -j 40 -p

Params: 
-s the absolute path of this snakemake script.
-j CPU core number
-p print the cmd used in every rule.
```

Wait until the workflow ends.



