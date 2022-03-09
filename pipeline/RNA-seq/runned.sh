#!/bin/bash

# ----- downloading public data -----
cd $workdir
download_raw_seq.py -p $PRJNAXXXXXX -f 'SRRxxxxxx1|SRRxxxxxx2|SRRxxxxxx3' -o 0_raw_data
cd ${workdir}/0_raw_data
bash ${PRJNAXXXXXX}.sh

# ----- run -----
cd $workdir
RNAseq_PE_step.sh $name 1,2,3,4 $threads $genomeVersion sample${i}_rep${j}.R1.fastq.gz sample${i}_rep${j}.R2.fastq.gz <ensGene|refGene> > bin/log/${name}.log 2>&1


# ----- extract mapping ratio -----

### content in file 'sample_list.txt' (seperated by tabs):
### $name sample${i}_rep${j}.R1.fastq.gz sample${i}_rep${j}.R2.fastq.gz

cd ${workdir}/bin
sed -i "s/ /\t/g" sample_list.txt
extract_log_fitCR.py sample_list.txt 'RNA-seq' data_information.csv
