# $ cat 225.bed |awk '{if($6=="+") print $1":"$3":"$6; else print $1":"$2":"$6}'|sort |uniq -c >225.uniq.pas
# shell too slow, using Python

# step1, get uniq pas number and freq
import re,time
# ## uniq
fIn=snakemake.input[0]
fOut=snakemake.output[0]

def now():
    return '['+ time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()) +']';

fr=open(fIn,'r')
i=0
DB={}
start=time.time();
while True:
    lineR=fr.readline();
    if not lineR:
        break;
    i+=1;

    if i%5e6==0:
        print( now() ,i,' line procced.')
    line=lineR.strip()
    arr=re.split("\t",line)
    #
    k=arr[0]+":"+arr[5]
    if arr[5]=='+':
        pos=arr[2]
    else:
        pos=arr[1]
    if k not in DB:
        DB[k]={};
    if pos not in DB[k]:
        DB[k][pos]=0
    DB[k][pos]+=1;
fr.close()
print('==========>line number: ', i)


# step2，save uniq PAS sites
# bed file column 5th:count, 2nd-3rd: the same position.
# chr1	564599	564599	chr1:564599:+	3.33339746426	+
fw=open(fOut,'w')
i=0
for chrStrand in DB:
    arr=re.split(':',chrStrand);
    #print(arr_chrStrand, len(DB[chrStrand])) #['chr10', '-'] 61299
    for pos in DB[chrStrand]:
        i+=1;
        count=DB[chrStrand][pos];
        arr2=[arr[0], pos,pos, arr[0]+":"+pos+':'+arr[1], str(count), arr[1]] 
        #['chr10', '101688', '101688', 'chr10:101688:-', '1', '-']
        fw.write('\t'.join(arr2) +'\n');
fw.close()
