#a06_addFreq_to_pas.py


fInBed=snakemake.input[0]
fInRefer=snakemake.input[1]
fOut=snakemake.output[0]


##########
## step1: load small bed of each cell, to a dict bedDB
##########
# key is chr1:+ 
# value is an array of dict {1:count1, 29:count2, ...}
import re,time

def now():
    return '['+ time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()) +'] ';


def readBed(fInBed):
    start=time.time();
    fr=open(fInBed,'r');
    i=0
    bedDB={}
    for lineR in fr.readlines():
        i+=1
        line=lineR.strip();
        arr=re.split('\t', line);
        #print(arr) #['chr10', '300469', '300566', 'ST-E00167:377:H52JTALXX:8:2202:22039:22194', '255', '+']
        #            ['chr10', '856735', '856853', 'ST-E00167:377:H52JTALXX:8:2216:8410:1608', '255', '-']
        #create (key)chrStrand in bedDB
        chrStrand=arr[0]+":"+arr[-1]
        if chrStrand not in bedDB:
            bedDB[chrStrand]={}
        #
        if arr[-1]=="+":
            pos =arr[2];
        else:
            pos=arr[1];
        #whether pos in chr1:+:{ this }
        if pos in bedDB[chrStrand]:
            bedDB[chrStrand][pos]+=1;
        else:
            bedDB[chrStrand][pos]=1
    fr.close();
    
    print(now(), 'readBed(): i=',i,'; ',' consume ',round(time.time()-start,2),'s; ', sep='')
    # readBed(): i=1223251;  2.57s; 20190928-222958
    return (i,bedDB);

lineNo, bedDB1 = readBed(fInBed);


# load referDB JSON
import json
with open(fInRefer,'r') as load_f:
    pasPostions = json.load(load_f)
    #print(type(pasPostions)) #<class 'dict'> 
#
print(now(), "===> step 1 done!")




##########
## step2: Based on small bedDB, refer to referDB, get counts for each row in referDB
##########
def savePASbyID(fOut,lineNo,bedDB1):
    start= time.time()
    i=0
    j=0 #counter
    fw=open(fOut,'w')
    for chrStrand in pasPostions:
        chr1,strand=re.split(':',chrStrand);
        for site in pasPostions[chrStrand]:
            i+=1;
            key=chr1+":"+site+":"+strand;
            count=0;
            #print('pasPostions[chrStrand][site]=',pasPostions[chrStrand][site])

            for pos in pasPostions[chrStrand][site]:
                pos=str(pos);
                if chrStrand in bedDB1:
                    if pos in bedDB1[chrStrand]:
                        count+=bedDB1[chrStrand][pos]
            fw.write('\t'.join([key, str(count)]) + '\n')
            j+=count;
    fw.close();

    print(now(), 'savePASbyID(): PAS site number i=',i, '; cover counts j=',j,'(',round(j/lineNo*100,3),'%); Consume time=',round(time.time()-start,2),'s; ', sep='')
savePASbyID(fOut,lineNo, bedDB1)


print(now(), "===> step 2 done!")

