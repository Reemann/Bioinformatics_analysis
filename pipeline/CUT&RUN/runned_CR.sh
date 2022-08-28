#!/bin/bash

# ----- downloading public data -----
cd $workdir
download_raw_seq.py -p $PRJNAXXXXXX -f 'SRRxxxxxx1|SRRxxxxxx2|SRRxxxxxx3' -o 0_raw_data
cd ${workdir}/0_raw_data
bash ${PRJNAXXXXXX}.sh


# ----- Pair-end -----
cd $workdir
CUT_RUN_TF_step.sh $name 1,2,3 $threads $genomeVersion sample${i}_rep${j}.R1.fq.gz sample${i}_rep${j}.R2.fq.gz true >> bin/log/${name}.log 2>&1


# ----- FRiP -----
cd $workdir

# private peaks
# This is generated defautly by CUT_RUN_TF_step.sh; It's not necessary to do it again!
# common_peak_FRiP.sh $name PE CR 3_peak/${name}_peaks.narrowPeak

# common peaks
path_to_common_peak_BED=path/to/common_peaks.bed
common_peak_FRiP.sh $name PE CR $path_to_common_peak_BED


# ----- BLAST low mapping samples -----
mapping_log_dir=path/to/MappingLogs
fastq_dir=path/to/fastqFiles # must only contains the fq to blast, all fq under this dir will be blast
blast_db=~/source/blast_db
# the order of mapping.logs must be consistent with the alphanumeric `ls` order of fastq files
blast_low_mapped_samples.py -d $fastq_dir -b $blast_db -n [$select_reads_number] -@ [$blast_threads] -m ${mapping_log_dir}/Mapping_${name1}.log,${mapping_log_dir}/Mapping_${name2}.log,${mapping_log_dir}/Mapping_${name3}.log


# ----- sampling down -----
# number of fragments should be modified
normalization=5000000
function processing {
    name=${1}
    # number of reps should be modified
    cat ../2_signal/${name}_rep[1-3]_fragments.bed > ../2_signal/${name}_merged_total_fragments.bed
    n=`cat ../2_signal/${name}_merged_total_fragments.bed | wc -l `
    ratio=`bc -l <<< "${normalization} / ${n}"`
    echo ${ratio}
    cat ../2_signal/${name}_merged_total_fragments.bed | awk -v ratio=${ratio} 'BEGIN{{srand(1006)}} {{if(rand()<ratio) print $0}}' > ../2_signal/${name}_total_fragments.bed

    macs2 callpeak -f BEDPE -t ../2_signal/${name}_total_fragments.bed -n ${name} -g 1.4e9 -q 0.01 --outdir ./ --keep-dup all 2>&1 >>/dev/null | tee ${name}_total_MACS.out

    fragment_length=`awk 'BEGIN{s=0} {s+=$3-$2} END{printf "%f", s/NR}' ../2_signal/${name}_total_fragments.bed`
    /mnt/Storage/home/wangyiman/bin/utilities/ShiftPairEnd.sh ../2_signal/${name}_total_fragments.bed ../2_signal/${fragment_length}
    n=`wc -l ../2_signal/${name}_total_fragments_shift.bed | cut -f 1 -d " "` && \
    c=`bc -l <<< "1000000 / $n"` && \
    genomeCoverageBed -bga -scale $c -i ../2_signal/${name}_total_fragments_shift.bed -g ~/source/bySpecies/danRer11_2/chromSizes/danRer11_2.chrom.sizes > ../2_signal/${name}_total_fragments_shift.bdg && \
    /mnt/Storage/home/wangyiman/bin/utilities/bdg2bw.sh ../2_signal/${name}_total_fragments_shift.bdg ~/source/bySpecies/danRer11_2/chromSizes/danRer11_2_main.chrom.sizes ${name} && \
    rm ../2_signal/${name}_total_fragments_shift.bed ../2_signal/${name}_total_fragments_shift.bdg
}

processing $sample



# ----- filter peaks -----
out_dir=${out_dir}
qvalue_cutoff=10
fold_cutoff=10

function filter {
    cd 3_peak
    sample=${1}
    awk -v fold_cutoff=${fold_cutoff} -v qvalue_cutoff=${qvalue_cutoff} 'NR>21{if($8>fold_cutoff && $9>qvalue_cutoff) printf "%s\t%d\t%d\t%s\t%d\t.\t%.5f\t%.5f\t%.5f\t%d\n", $1, $2-1, $3, $10, $9*10, $8, $7, $9, $5-$2}' ${sample}_peaks.xls > ${out_dir}/${sample}_filtered_peaks.narrowPeak
    awk -v fold_cutoff=${fold_cutoff} -v qvalue_cutoff=${qvalue_cutoff} 'NR>21{if($8>fold_cutoff && $9>qvalue_cutoff) printf "%s\t%d\t%d\t%s\t%.5f\n", $1, $5-1, $5, $10, $9}' ${sample}_peaks.xls > ${out_dir}/${sample}_filtered_summits.bed
    cd ../
}

cd $workdir
filter sample

