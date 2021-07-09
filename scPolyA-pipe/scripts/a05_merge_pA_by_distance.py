################
# Aim: recode, create APA index file
# merge APA sites within 24/2=12nt, get feature row name.
# if cluster longer than 24nt, then 12nt before and after leaves, others recluster.
################

fIn=snakemake.input[0];
fOut=snakemake.output[0];

import sys
fLog=open(snakemake.log[0], "w")
sys.stderr = sys.stdout = fLog

###############
## step1: Read to dict, if count>1
###############
#file name
import re,time
#chr10	159177	159177	chr10:159177:-	20	-
fr=open(fIn,'r')
start=time.time();

def now():
    return '['+ time.strftime("%Y/%m/%d %H:%M:%S", time.localtime()) +'] ';


i=0;
j=0;
pasDB={} #value is array chr1:+ : [1,2,3] position array. May change later, keep position within cluster
pasDict={} #value is dict(position: counts)  chr1:+ : {1: 21, 3:45, 10:2, } pos and counts, not change
for lineR in fr.readlines():
    i+=1
    line=lineR.strip();
    arr=re.split('\t',line) #['chr10', '101688', '101688', 'chr10:101688:-', '1', '-']
    #filter: by count. May keep??//todo
    if int(arr[4])<=1:
        continue;
    j+=1;
    k=arr[0]+':'+arr[-1]
    if k not in pasDB:
        pasDB[k]=[]
        pasDict[k]={}
    pasDB[k].append(int(arr[1]));
    pasDict[k][arr[1]]=int( arr[4] )
fr.close();

# log info
print(now(), "length of pasDB and pasDict", len(pasDB), len(pasDict));
for k in pasDB.keys():
    print(k, len(pasDB[k]) )
print(now(), 'Total PAS number i=', i, '; PAS number who has enough support number j=', j,'; Ratio j/i=', round(j/i*100),'%', sep='') #2732224

# Total PAS number i=2732224; PAS number who has enough support number j=1574173; Ratio j/i=58%
print(now(), "step 1 done!")






###############
## step2: within same chr and strand, sort, cluster by distance.
###############

def getClusterFromPasDB(pasDB, distance=24):
    pasCluster={} #total cluster, save cluster
    i=0
    start=time.time();
    for chrStrand in pasDB: #chr1:+  : [ 1,2,3,100,...]
        pasDB[chrStrand].sort();# order
        pasCluster[chrStrand]=[]; #accept rs after cluster  chr1:+  : [ [1,2,3],  [100,113,120], ]
        tmp=[]; # current cluster pos. In the end, add whole cluster to the array.
        lastPos=-1;
        for pos in pasDB[chrStrand]:
            i+=1;
            if lastPos==-1: #when init
                tmp.append(pos); 
                lastPos=pos;
                continue;
            # if within <=distance, go on merge
            if pos-lastPos<=distance:
                tmp.append(pos); 
                lastPos=pos;
            else: # if > distance, then gave tmp to array cluster. Clean tmp for new commer next.
                pasCluster[chrStrand].append(tmp);
                tmp=[pos];
                lastPos=pos;
        if len(tmp)>0:#fix: the remaining of the chr is a cluster
            pasCluster[chrStrand].append(tmp);
    #
    print(now(), "getClusterFromPasDB(): Total row number i=",i, '; distance=',distance,'; ',' consume time:',round(time.time()-start,2),"s", sep='')
    #==end== Total row number i=1574173; distance=24; consume time0.84;
    return pasCluster;
#
distance=24 #filter by distance
pasCluster1=getClusterFromPasDB(pasDB, distance)
#check()
print(now(), "step 2 done!")





###############
## step3 [change some variables, so start from step2 if redo] 
###############
# Parametr: distance Or distance/2, for long Cluster?

# For cluster span more than 24nt, get the Pos with most counts
# Check pasDict, based on Pos array and chrStrand, get the most reads, return the Pos
def getMaxPos(win, chrStrand):
    maxPos=win[0]
    maxCount=pasDict[chrStrand][str(maxPos)]
    #
    for pos in win:
        count=pasDict[chrStrand][str(pos)];
        #print(i, count)#133781868 804
        if count>maxCount:
            maxCount=count;
            maxPos=pos;
    #print(maxPos, maxCount)
    return maxPos
#test
#win=[133781863, 133781864, 133781865, 133781866, 133781867, 133781868, 133781869, 133781870, 133781884, 133781887, 133781888]
#getMaxPos(win,'chr10:-') #133781868


# keep pos within distance, others put back to pasDB2, and retun.(in order to the 2nd cluster)
def getPosFarFromMaxPos(pasCluster, distance):
    i=0
    j=0
    start=time.time()
    pasDB2={} #value is arr chr1:+ : [1,2,3] Sites whose distance to top counts too far

    for chrStrand in pasCluster:
        #chrStrand='chr10:-';
        pasDB2[chrStrand]=[]
        for index in range( len(pasCluster[chrStrand]) ):
            win=pasCluster[chrStrand][index]
            i+=1;
            #recluster span longer than 24nt
            if max(win)-min(win)>distance:
                inWin=[] #accept within, at last overwrite win
                maxPos=getMaxPos(win,chrStrand); # find the pos of maxCount, then find distance range
                for pos in win:
                    if (maxPos-distance/2)<=pos<=(maxPos+distance/2): ##>>>>>>>>>>> critical
                        inWin.append(pos);
                    else:
                        pasDB2[chrStrand].append(pos); # if not in range, give back to DB, then recluster
                #print('chr10:'+str( min(win) )+'-'+str(max(win)),   '\t',length,'nt', win )
                # chr10:320394-320426 	 32 nt [320394, 320395, 320396, 320397, 320398, 320399, 320400, 320401, 320402, 320403, 320404, 320405, 320406, 320407, 320410, 320411, 320412, 320414, 320422, 320423, 320426]
                # chr10:856724-856762 	 38 nt [856724, 856731, 856732, 856733, 856734, 856735, 856736, 856737, 856738, 856739, 856740, 856741, 856742, 856744, 856745, 856761, 856762]
                j+=1;
                pasCluster[chrStrand][index]=inWin; #win=inWin;
    #
    print(now(), 'getPosFarFromMaxPos():Total cluster number i=',i, '; Too long cluster number j=',j, '; Ratio j/i=', round(j/i*100,2),'% ','; Time consuming:',round(time.time()-start,2),'s; ', sep='')
    # i= 25600 ; j= 1001 ; j/i= 4 % cluster 24nt
    # i= 16237 ; j= 439 ; j/i= 2.7 % 
    # i= 16400 ; j= 351 ; j/i= 2 % #distance=20
    # i= 17144 ; j= 176 ; j/i= 1 % #distance=10
    # i= 18510 ; j= 79 ; j/i= 0.43 % #distance=5
    #getPosFarFromMaxPos():Total cluster number i=708828; Too long cluster number j=21201; Ratio j/i=2.99% ; Time consuming:0.61s;
    return pasDB2;
#
pasDB2=getPosFarFromMaxPos(pasCluster1, distance)

print(now(), "step 3 done!")





###############
## step4 begin loop
###############
# get the number of pos in the value array of current cluster.
# chr1:+ :[1,2,300, ...]
def getPosNumber(pasDB2):
    number=0
    for chrStrand in pasDB2:
        arr=pasDB2[chrStrand];
        number+=len(arr)
    #print(number)
    return number;
#
#print('pasDB: ',getPosNumber(pasDB) )
#print('pasDB2: ', getPosNumber(pasDB2))

def getClusterNumber(pasCluster):
    number=0
    for chrStrand in pasCluster:
        arr=pasCluster[chrStrand];
        number+=len(arr)
    #print(number)
    return number;


roundNum=0
while True:
    number=getPosNumber(pasDB2);
    roundNum+=1;
    print(now(), '='*30 , 'roundNum=',roundNum, '; number left ', number);
    
    if number==0 or roundNum>100:
        break;
    
    #go to step2, using pasDB2 recluster
    pasCluster2=getClusterFromPasDB(pasDB2, distance)
    
    #go to step3, keep sites widin distance nt of maxPos, others back to pasDB2.
    pasDB2=getPosFarFromMaxPos(pasCluster2, distance)
    
    # merge pasCluster2 to pasCluster1
    i=0
    for chrStrand in pasCluster2:
        for cluster in pasCluster2[chrStrand]:
            i+=1
            pasCluster1[chrStrand].append(cluster);
            #print(chrStrand,cluster)
    print(now(), i, 'cluster combined to pasCluster1, Total cluster number: %d\n' % getClusterNumber(pasCluster1) )

print(now(), "step 4 done!")





###############
## step5: check total number of clusters
###############
w=0
for chrStrand in pasCluster1:
    i=0;
    for win in pasCluster1[chrStrand]:
        i+=1;
    print(chrStrand, i)
    w+=i;
print(now(), "[Check] Total cluster number: ", w) #742288 clusters


# check how many site in the cluster
#chr1:+  : [ [1,2], [300,314],... ]
def getPosNumberFrom(pasCluster1):
    i=0
    for chrStrand in pasCluster1:
        for c in pasCluster1[chrStrand]:
            if len(c)==0:
                print(now(), "Empty window: ", c )
            i+=len(c)
    return i;
    #print(i)
print(now(), 'pasCluster1:', getPosNumberFrom(pasCluster1) ) #1569257

def check():
    n1=getPosNumber(pasDB2);
    n2=getPosNumberFrom(pasCluster1);
    print(now(), 'Check', '>'*10, "pasDB2:%d, pasCluster1:%d; Total:%d" %(n1, n2, n1+n2))
check()


# step2.5 check result
# print() cluster result, who span over 24 nt
def getLongClusterNumber(pasCluster, distance):
    i=0
    j=0
    pasDB2={} #value like   arr chr1:+ : [1,2,3]
    tmpDB={};

    for chrStrand in pasCluster:
        #chrStrand='chr10:-';
        tmpDB[chrStrand]=[];
        for win in pasCluster[chrStrand]:
            i+=1;
            length=max(win)-min(win);
            if length>distance:
                if j<5: #print some examples
                    print('Warning: too long window: chr10:'+str( min(win) )+'-'+str(max(win)),   '\t',length,'nt', win )
                j+=1;
    #
    print(now(), 'cluster total number i=',i, '; cluster which longer than ',distance,'nt j=',j, '; Long/total cluster ratio j/i=', round(j/i*100,2),'%', sep='')
    # i= 25600 ; j= 1001 ; j/i= 4 % region, cluster span > 24nt
    # i= 16237 ; j= 439 ; j/i= 2.7 % filtr out  site <=1
    # i= 16400 ; j= 351 ; j/i= 2 % #distance=20
    # i= 17144 ; j= 176 ; j/i= 1 % #distance=10
    # i= 18510 ; j= 79 ; j/i= 0.43 % #distance=5
    return j;
#
getLongClusterNumber(pasCluster1, distance)

print(now(), "step 5 done!")





###############
## step6: get pos with most Counts within each cluster, save to another array: pasPostions
###############
# key is chr1:+ ， value is a dit {1:[1,24],  121:[100,121],  300:[300], ...} 
#                         top count's pos:[arr of all pos in the cluster]
## for debug
def print2(*args):
    if printFlag==False: #False output nothing, True normal
        return;
    print(args)
#test
printFlag=False;

start=time.time();
pasPostions={}
i=0
for chrStrand in pasCluster1:
    pasPostions[chrStrand]={};
    for cluster in pasCluster1[chrStrand]:
        i+=1;
        maxPos=getMaxPos(cluster,chrStrand);
        
        pasPostions[chrStrand][maxPos]=cluster;# pos of maxCount as key, value is all the pos in the cluster
    #print(chrStrand, i)
print(now(), i,' clusters; Time consuming:', round(time.time()-start,2),'s' , sep='') #742288 clusters
print(now(), "step 6 done!")




###############
## step7 save the array pasPostions to file
###############
import json

#json to file
with open(fOut, "w") as fw:
    json.dump(pasPostions, fw)
# 26M Sep 27 17:03 pasPostions_siteAndIncludes_PY.json
# less pasPostions_siteAndIncludes_PY.json

print(now(), "step 7 done!")








###############
## step8 output pasPostions to bed format, 
###############
# to compare with polyA_DB2/3, or filter and validation
# chr10	159177	159177	chr10:159177:-	20	-

fOutBed=snakemake.output[1]

fwBed=open(fOutBed,"w")
start=time.time();
i=0
flag=True;
for chrStrand in pasPostions:
    #pasPostions str
    chrN,strand=re.split(':', chrStrand);
    # key is chr1:+ ， value is a dict {1:[1,24],  121:[100,121],  300:[300], ...} 
    # Pos of maxCount:[ all the pos in the cluster]
    for site in pasPostions[chrStrand]:
        i+=1;
        cluster=pasPostions[chrStrand][site];
        count=0;
        for pos in cluster:
            count+=pasDict[chrStrand][str(pos)];
        # to bed format
        site=str(site)
        arr=[chrN, site, site, chrN+":"+site+":"+strand, str(count), strand];
        fwBed.write('\t'.join(arr)+'\n');
    #print(chrStrand, i)
print(now(), i,' clusters; Time consuming:', round(time.time()-start,2),'s' , sep='')
fwBed.close()

print(now(), "step 8 done!")

# 35M Sep 27 17:33 pasPostions_siteAndIncludes_PY.bed
# $ head pasPostions_siteAndIncludes_PY.bed
#chr10	157706	157706	chr10:157706:-	5	-
#chr10	159177	159177	chr10:159177:-	20	-
#chr10	167933	167933	chr10:167933:-	3	-
#chr10	174098	174098	chr10:174098:-	142	-
#chr10	185668	185668	chr10:185668:-	13	-

# 742288clusters; Time consuming:1.75s




print(now(), "we get files: ", fOut, fOutBed)

fLog.close()

