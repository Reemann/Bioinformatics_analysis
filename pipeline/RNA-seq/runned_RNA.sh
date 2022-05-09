#!/bin/bash

# ----- downloading public data -----
cd $workdir
download_raw_seq.py -p $PRJNAXXXXXX -f 'SRRxxxxxx1|SRRxxxxxx2|SRRxxxxxx3' -o 0_raw_data
cd ${workdir}/0_raw_data
bash ${PRJNAXXXXXX}.sh


# ----- run -----
cd $workdir
# normal
RNAseq_PE_step.sh $name 1,2,3,4 $threads $genomeVersion sample${i}_rep${j}.R1.fastq.gz sample${i}_rep${j}.R2.fastq.gz <ensGene|refGene> > bin/log/${name}.log 2>&1
# calculate expression level by Exon
RNAseq_PE_step_ExonUniq.sh $name 1,2 $threads $genomeVersion sample${i}_rep${j}.R1.fastq.gz sample${i}_rep${j}.R2.fastq.gz <ensGene|refGene> > bin/log/${name}.log 2>&1
# if stringtie ended with segment fault and ended at a certain bundle ("chrXX:XXXXXXX-XXXXXXX") each time :
function rerun_segment_fault_samples{
    name=$1
    segment_fault_site=$2
    
    cd ~/maternal_loading/2.public_data/1_mapping
    samtools sort -@ 20 -o ${name}_sorted.bam ${name}.bam && samtools index -@ 20 ${name}_sorted.bam
    samtools view ${name}_sorted.bam $site > ${name}_sorted_segmentFaultRead.sam
    cut -f 1 ${name}_sorted_segmentFaultRead.sam | sort | uniq > ${name}_sorted_segmentFaultRead.list.txt
    java -jar ~/software/picard.jar FilterSamReads I=${name}_sorted.bam O=${name}_sorted_filtered.bam READ_LIST_FILE=${name}_sorted_segmentFaultRead.list.txt FILTER=excludeReadList

    cd ~/maternal_loading/2.public_data/2_expression_value_ExonUniq
    stringtie ../1_mapping/${name}_sorted_filtered.bam -p 30 -G /mnt/Storage/home/wangyiman/source/bySpecies/danRer11_2/ensGene/danRer11_2.ensGene.ExonUniq.gtf -l ${name} -A RNA_seq_ventricle_rep1_ensGene_coverage.ExonUniq.txt -o RNA_seq_${name}_ensGene_coverage.ExonUniq.gtf -e -B -v > ../bin/log/${name}_rerun.log 2>&1
    
    tail ../bin/log/${name}_rerun.log
}
rerun_segment_fault_samples $name "chrXX:XXXXXXX-XXXXXXX"


# ----- extract mapping ratio -----

### content in file 'sample_list.txt' (seperated by tabs):
### $name sample${i}_rep${j}.R1.fastq.gz sample${i}_rep${j}.R2.fastq.gz

cd ${workdir}/bin
sed -i "s/ /\t/g" sample_list.txt
extract_log_fitCR.py sample_list.txt 'RNA-seq' data_information.csv


# ----- view reads in IGV -----
cd ${workdir}/2_expression_value
samtools sort -@ $((${processer}-1)) -o ${name}_sorted.bam ${name}.bam && samtools index -@ $((${processer}-1)) ${name}_sorted.bam
