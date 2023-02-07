pyPATH=/mnt/Storage/home/yanghui/scripts/python

bedF=${1}
cgF=${2}
outF=${3}
threads=${4} # 6

# step1. split into multiple(=threads) small files with 
cut -f 1-3 ${bedF} | sort -u > tmp.bed # 
n_row=`cat tmp.bed | wc -l`
para_l=`bc <<< ${n_row}/${threads}+1`
split -l ${para_l} tmp.bed -d -a 3 ${outF}_subfile_

# step2. calculate average DNA methylation for each small file
for file in `ls ${outF}_subfile_*`
do
	grep '^chr' ${file} | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3}' | sort -u > ${file}.tmp && \
	intersectBed -wao -a ${file}.tmp -b ${cgF} > ${file}.all.G && \
	python ${pyPATH}/methyl_processing.py region_methylation ${file}.all.G ${file}.ave &
done # for file end
wait;
# cat ${outF}_subfile_*.ave | sort -S100G -k1,1 -k2,2n --parallel=6 > ${outF}
cat ${outF}_subfile_*.ave > ${outF}.tmp
rm ${outF}_subfile_*

# order according to primary input files
cut -f 1-3 ${bedF} | intersectBed -wao -a - -b ${outF}.tmp -f 1.00 -r | cut -f 1-3,7 > ${outF}
rm tmp.bed ${outF}.tmp
