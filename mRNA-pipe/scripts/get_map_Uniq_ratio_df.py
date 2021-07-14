fInArr=snakemake.input;
fOut=snakemake.output[0]
cids=snakemake.params['cids']
log_path=snakemake.params["log_path"]

import re

def getInfo(cid):
    fr=open(log_path+'/'+cid+"_Log.final.out")
    out=[cid]
    for lineR in fr.readlines():
        line=lineR.strip()
        #print(line)
        if line.startswith("Uniquely") or line.startswith("Number of input reads"):
            arr=re.split('\|\t', line)
            out.append(arr[1])
    fr.close()
    # cid total uniqReads uniqRatio
    return(out)


fw=open(fOut, 'w')
fw.write("cid\ttotal\tuniqReads\tuniqRatio\n")
for cid in cids:
  print('Record uniq mapped reads of cid = ', cid)
  arr=getInfo(cid)
  fw.write('\t'.join(arr)+"\n")
fw.close()
