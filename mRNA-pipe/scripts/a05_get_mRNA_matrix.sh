# must input 2 paras: idListFile outputFile [dir]
# aim: according id in idListFile, combine 2nd col of each file in dir, get matrix, output to outputFile
# note: user must insure rownames are the same. The script only combine the 2nd col, not order
# v3.0 for APA
# v3.1 for mRNA

## settings
idListFile=$1; # the first para is: id list file
idListFile=${idListFile:?"Must input: id list file path"}

outputFile=$2; #output filename
outputFile=${outputFile:?"Must input output filename"}

dir=$3; #small data file dir, default ./
dir=${dir:='.'}


## current temp dir name
tmp="tmp_"`date |sha256sum|base64|head -c 20`
tmp_test2="test2_"`date |sha256sum|base64|head -c 20`

##step1 add 1st column to target
id=`head -n 1 ${idListFile}`
awk '{print $1}' ${dir}/${id}.bam.count > $outputFile;

##step2 combine all 2nd columns to target
i=0;
cat ${idListFile} | while read id; do ((i++)); echo ${i} $id;
	awk '{print $2}' ${dir}/${id}.bam.count > $tmp;
	paste $outputFile $tmp >$tmp_test2;
	mv $tmp_test2 $outputFile;
done;

##step3 add the first row to target
line1=`cat ${idListFile} | xargs`;
line1="gene "${line1}
## echo ${line1};
sed -i "1i ${line1}" $outputFile

##step3.1 remove last 5 lines
grep -v "^_" ${outputFile} >$tmp 
mv $tmp ${outputFile}

#step4 clean
rm $tmp
echo "done~ result file: "${outputFile}
