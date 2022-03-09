#!/bin/bash

# ----- downloading public data -----
cd $workdir
download_raw_seq.py -p $PRJNAXXXXXX -f 'SRRxxxxxx1|SRRxxxxxx2|SRRxxxxxx3' -o 0_raw_data
cd ${workdir}/0_raw_data
bash ${PRJNAXXXXXX}.sh


# ----- Single-end narrow peaks -----
cd $workdir
ChIP_SE_nucleosome_narrow_step.sh $name 1,2,3 $threads $genomeVersion sample${i}_rep${j}fastq.gz true >> bin/log/${name}.log 2>&1 &


# ----- Pair-end narrow peaks -----
cd $workdir
ChIP_PE_nucleosome_narrow_step.sh $name 1,2,3 $threads $genomeVersion sample${i}_rep${j}.R1.fastq.gz sample${i}_rep${j}.R2.fastq.gz true >> bin/log/${name}.log 2>&1 &


