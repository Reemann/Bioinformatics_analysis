#!/bin/bash

cd $workdir
RNAseq_PE_step.sh $name 1,2,3,4 $threads $genomeVersion sample${i}_rep${j}.R1.fastq.gz sample${i}_rep${j}.R2.fastq.gz <ensGene|refGene> > bin/log/${name}.log 2>&1

### content in file 'sample_list.txt' (seperated by tabs):
### $name sample${i}_rep${j}.R1.fastq.gz sample${i}_rep${j}.R2.fastq.gz

cd ${workdir}/bin
sed -i "s/ /\t/g" sample_list.txt
extract_log_fitCR.py sample_list.txt 'RNA-seq' data_information.csv
