pyPATH=/mnt/Storage/home/yanghui/scripts/python

bedF=${1}
bwF=${2}
n=${3}
outF=${4}

grep '^chr' ${bedF} | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3}' > ${outF}.bed.tmp
python ${pyPATH}/getBigWigValue.py ${outF}_tmp ${outF}.bed.tmp ${n} ${bwF} # ${outF}_tmp_siteprof${n}[temp]
gunzip -f ${outF}_tmp_siteprof1.gz
cut -f 4- ${outF}_tmp_siteprof1 > ${outF}
rm ${outF}_tmp_siteprof1 ${outF}.bed.tmp