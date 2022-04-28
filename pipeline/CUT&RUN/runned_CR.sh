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
fastq_dir=path/to/fastqFiles
blast_db=~/source/blast_db
# the order of mapping.logs must be consistent with the alphanumeric `ls` order of fastq files
blast_low_mapped_samples.py -d $fastq_dir -b $blast_db -n [$select_reads_number] -@ [$blast_threads] -m ${mapping_log_dir}/Mapping_${name1}.log,${mapping_log_dir}/Mapping_${name2}.log,${mapping_log_dir}/Mapping_${name3}.log
