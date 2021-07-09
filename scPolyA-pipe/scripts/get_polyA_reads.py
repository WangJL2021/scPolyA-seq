# Aim: get polyA reads:
# keep 10A or 6A$, remove too short reads(<=10)

import re, time;

def getPolyARead(fileNameArr):
    # read file
    fr=open(fileNameArr[0], "r")
    # write to files:  keep(10A /end 6A), other
    fwKeep=open(fileNameArr[1], 'w')
    fwOther=open(fileNameArr[2], 'w')

    i=0
    fastq=[] #fq file 4 row
    tag="" #tag for each 4 row
    start=time.time()
    for lineR in fr.readlines():
        i+=1
        line=lineR.strip()
        fastq.append(lineR)

        #give tag by scan sequence (2nd row)
        if i%4==2:
            haveA10=re.search("A{10,}", line)
            haveA6=re.search("A{6,}$", line)
            #
            if haveA10:
                tag="10A"
            elif haveA6:
                tag="Tail6A"
            else:
                tag="Other"

        # at row 4, output by tag
        if i%4==0:
            # save to mem by tag at row 1,3
            if tag=="10A" and (haveA10.span()[0]>10):
                #print(haveA10.span()[0],"\n")
                fastq[1]=fastq[1].strip()[0: haveA10.span()[0] ]+"\n"
                fastq[3]=fastq[3].strip()[0: haveA10.span()[0] ]+"\n"
                for item in fastq:
                    fwKeep.write(item)
            elif tag=="Tail6A" and (haveA6.span()[0]>10):
                #print(haveA6.span()[0],"\n")
                fastq[1]=fastq[1].strip()[0: haveA6.span()[0] ]+"\n"
                fastq[3]=fastq[3].strip()[0: haveA6.span()[0] ]+"\n"
                for item in fastq:
                    fwKeep.write(item)
            else:
                # too short, or no 10A or 6A
                for item in fastq:
                    fwOther.write(item)
             #clen
            fastq=[]
            
    #close fiels.
    fr.close()
    fwKeep.close()
    fwOther.close()
    return 0;

#
fileNameArr=[snakemake.input[0], snakemake.output[0], snakemake.output[1]]
getPolyARead(fileNameArr)

