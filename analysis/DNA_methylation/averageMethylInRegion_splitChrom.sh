#!/bin/bash

pyPATH=/mnt/Storage/home/yanghui/scripts/python

bedF=${1}
cgF=${2}
outF=${3}
genome=${4}

case ${genome} in
    "hg38")
        chr_max=22
        ;;
    "mm10")
        chr_max=19
        ;;
    *)
        echo "Error!! Please provide a genome version!!"
        exit
esac

bedF_base1=$(basename ${bedF})
bedF_base=${bedF_base1%%.bed}
cgF_base1=$(basename ${cgF})
cgF_base=${cgF_base1%%.bed}

cgF_chrs_dir=${outF%%.*}_tmp
mkdir -p ${cgF_chrs_dir}

# step1. split by chromosome
### 1 - 5
for chrN in 1 2 3 4 5 X Y;do
    echo $chrN
    [[ ! -e ${bedF_base}.chr${chrN}.bed ]] && cut -f 1,2,3 ${bedF} | grep ^chr${chrN}$'\t' | sort -k1,1 -k2,2n > ${bedF_base}.chr${chrN}.bed
    [[ ! -e ${cgF_chrs_dir}/${cgF_base}.chr${chrN}.bed ]] && grep ^chr${chrN}$'\t' ${cgF} | sort -k1,1 -k2,2n > ${cgF_chrs_dir}/${cgF_base}.chr${chrN}.bed
    intersectBed -wa -wb -a ${bedF_base}.chr${chrN}.bed -b ${cgF_chrs_dir}/${cgF_base}.chr${chrN}.bed > ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.bed && \
    python ${pyPATH}/methyl_processing.py region_methylation ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.bed ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.ave.bed && \
    intersectBed -wao -a ${bedF_base}.chr${chrN}.bed -b ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.ave.bed -f 1.00 -r | cut -f 1-3,7 > ${cgF_chrs_dir}/${outF}.chr${chrN} &
done
wait

### 6 - 10
for chrN in 6 7 8 9 10;do
    echo $chrN
    [[ ! -e ${bedF_base}.chr${chrN}.bed ]] && cut -f 1,2,3 ${bedF} | grep ^chr${chrN}$'\t' | sort -k1,1 -k2,2n > ${bedF_base}.chr${chrN}.bed
    [[ ! -e ${cgF_chrs_dir}/${cgF_base}.chr${chrN}.bed ]] && grep ^chr${chrN}$'\t' ${cgF} | sort -k1,1 -k2,2n > ${cgF_chrs_dir}/${cgF_base}.chr${chrN}.bed
    intersectBed -wa -wb -a ${bedF_base}.chr${chrN}.bed -b ${cgF_chrs_dir}/${cgF_base}.chr${chrN}.bed > ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.bed && \
    python ${pyPATH}/methyl_processing.py region_methylation ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.bed ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.ave.bed && \
    intersectBed -wao -a ${bedF_base}.chr${chrN}.bed -b ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.ave.bed -f 1.00 -r | cut -f 1-3,7 > ${cgF_chrs_dir}/${outF}.chr${chrN} &
done
wait

### 11 - 15
for chrN in 11 12 13 14 15;do
    echo $chrN
    [[ ! -e ${bedF_base}.chr${chrN}.bed ]] && cut -f 1,2,3 ${bedF} | grep ^chr${chrN}$'\t' | sort -k1,1 -k2,2n > ${bedF_base}.chr${chrN}.bed
    [[ ! -e ${cgF_chrs_dir}/${cgF_base}.chr${chrN}.bed ]] && grep ^chr${chrN}$'\t' ${cgF} | sort -k1,1 -k2,2n > ${cgF_chrs_dir}/${cgF_base}.chr${chrN}.bed
    intersectBed -wa -wb -a ${bedF_base}.chr${chrN}.bed -b ${cgF_chrs_dir}/${cgF_base}.chr${chrN}.bed > ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.bed && \
    python ${pyPATH}/methyl_processing.py region_methylation ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.bed ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.ave.bed && \
    intersectBed -wao -a ${bedF_base}.chr${chrN}.bed -b ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.ave.bed -f 1.00 -r | cut -f 1-3,7 > ${cgF_chrs_dir}/${outF}.chr${chrN} &
done
wait

### 16 - max
for chrN in $(seq 16 1 ${chr_max});do
    echo $chrN
    [[ ! -e ${bedF_base}.chr${chrN}.bed ]] && cut -f 1,2,3 ${bedF} | grep ^chr${chrN}$'\t' | sort -k1,1 -k2,2n > ${bedF_base}.chr${chrN}.bed
    [[ ! -e ${cgF_chrs_dir}/${cgF_base}.chr${chrN}.bed ]] && grep ^chr${chrN}$'\t' ${cgF} | sort -k1,1 -k2,2n > ${cgF_chrs_dir}/${cgF_base}.chr${chrN}.bed
    intersectBed -wa -wb -a ${bedF_base}.chr${chrN}.bed -b ${cgF_chrs_dir}/${cgF_base}.chr${chrN}.bed > ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.bed && \
    python ${pyPATH}/methyl_processing.py region_methylation ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.bed ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.ave.bed && \
    intersectBed -wao -a ${bedF_base}.chr${chrN}.bed -b ${cgF_chrs_dir}/${bedF_base}__${cgF_base}.chr${chrN}.ave.bed -f 1.00 -r | cut -f 1-3,7 > ${cgF_chrs_dir}/${outF}.chr${chrN} &
done
wait

cat ${cgF_chrs_dir}/${outF}.chr* | sort -k1,1 -k2,2n > ${outF}

