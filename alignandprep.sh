#!/bin/bash
begin=$(date +"%s")

cd /wsu/home/fh/fh85/fh8591/piquelab/cindy/STARR/
#EXP=Bp1R14T6
#FASTQ=/nfs/rprscratch/Nextseq/170328_NS500258_0127_AHJNFWAFXX/fastqTrim3
#IDENTITY=032917
#SAMPLE=S1

mkdir -p ${IDENTITY}
mkdir -p ${IDENTITY}/working

cd ${IDENTITY}

less ${FASTQ}/${EXP}_${SAMPLE}_R2_001.fastq.gz | 
awk '{ 
	getline; print $0,"\n"" "; getline; print $0; getline; print $0}' | 
gzip > ${EXP}.I1singleline.fastq.gz 

#label read fastq
paste <(zcat ${FASTQ}/${EXP}_${SAMPLE}_R1_001.fastq.gz) <(zcat ${EXP}.I1singleline.fastq.gz) | awk '{
  if (NR%4 == 1)
  print $1"_"$3,$2;
else
  print $1;}' | gzip > ${EXP}_labeled.fastq.gz 

#check if there are short UMIs (<10bp, ususally see some small proportion of 8 and 9bp)
less ${FASTQ}/${EXP}_${SAMPLE}_R2_001.fastq.gz | awk '{if (NR%4==2) print length}' | sort | uniq -c > lengthUMI.${EXP}.txt &

#somtimes labeling stops before finished, this way can double check if the reads listed by hisat aren't right
less ${EXP}_labeled.fastq.gz | wc -l > length_labeled.${EXP}.txt & #1946453/4 fixed= /4

module load hisat2/2.0.4

hisat2 -x /wsu/home/groups/piquelab/fungei/atac/hisat.ref/grch37_snp/genome_snp --add-chrname -p 12 -U ${EXP}_labeled.fastq.gz --rg-id ${EXP} --rg SM:${EXP} | samtools view -@ 12 -b1 - > ${EXP}.bam 

samtools view -h ${EXP}.bam | 
  awk '{if (FNR<=87) print $0; else if($1 ~ /[.*_][ACT][CGT][ACG][AGT][ACT][CGT][ACG][AGT][ACT][CT]$/) print $0; }' | 
  samtools view -@ 12 -Sb1 > ${EXP}.f.bam 

samtools sort ${EXP}.f.bam -o ${EXP}.s.bam 
samtools index ${EXP}.s.bam 

module load gnu-4.7.2/python-2.7
export PYTHONPATH=/wsu/home/groups/piquelab/lib/python2.7/site-packages:$PYTHONPATH
export PYTHONPATH=/wsu/home/fh/fh85/fh8591/.local/lib/python2.7/site-packages:$PYTHONPATH

/wsu/home/fh/fh85/fh8591/.local/bin/umi_tools dedup -I ${EXP}.s.bam --mapping-quality=20 -S dedup.${EXP}.q20.bam -L dedup.${EXP}.log

#/wsu/home/fh/fh85/fh8591/.local/bin/umi_tools dedup -I ${EXP}.s.bam --output-stats=dedup.${EXP} -S dedup.${EXP}.bam 

#samtools view -h dedup.${EXP}.q20.bam | awk '{if (FNR<=87) print $0; else if($1 ~ /[.*_][ACT][CGT][ACG][AGT][ACT][CGT][ACG][AGT][ACT][CT]/) print $0;}' > dedup.${EXP}.q20.sam 
#samtools view -@ 12 -Sb1 dedup.${EXP}.q20.sam > dedup.${EXP}.q20_UMIF.bam 
samtools sort dedup.${EXP}.q20.bam -o dedup.${EXP}.q20.s.bam 
samtools index dedup.${EXP}.q20.s.bam 
#samtools view dedup.${EXP}.q20.bam  | wc -l > dedup.${EXP}_q20.counts.txt &
samtools view dedup.${EXP}.q20.bam  | wc -l > dedup.${EXP}_q20.counts.txt &
#call variants
samtools mpileup -f /wsu/home/groups/piquelab/data/RefGenome/hg19.fa dedup.${EXP}.q20.s.bam -l ~/piquelab/cindy/bitStarr/old/mpraSiteSelection/mpraSNPs_090215.reformated.sort.bed -t DP4 -g -d 1000000 > dedup.${EXP}.q20.pileup.bcf 
bcftools call -m -Oz dedup.${EXP}.q20.pileup.bcf > dedup.${EXP}.q20.call.bcf 
bcftools index dedup.${EXP}.q20.call.bcf 
bcftools query -f  '%CHROM\t%POS\t%REF\t%ALT{0}[\t%SAMPLE\t%DP4]\n' dedup.${EXP}.q20.call.bcf | awk -v OFS='\t' '{gsub(",","\t");print}' > working/dedup.${EXP}.q20.query.txt 



termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed for Script Execution. #Completed: ${REF} ${NAME}"