#!/bin/bash


######################
# Begin work section #
######################


## This script includes the following jobs:
# 1.  adapter trimming from raw data
# 2.  barcode extraction from R2
# 3.  (optional) distribution of clean read length
# 4.  rRNA mapping and the removal of mapped reads
# 5.  small RNA mapping and the removal of mapped reads
# 6.  genome mapping
# 7.  bam indexing and statistics of flags from the raw alignment file
# 8.  Remove unmapped, mate unmapped, not primary alignment, reads failing platform
# 9.  bam indexing and statistics of flags from the refined alignment file (step 8)
# 10. Remove PCR duplicates
# 11. bam indexing and statistics of flags from the dedup alignment file (step 10)
# 12. mutation calling using pileup2var


echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)


## main setting
# sample name: e.g., col_FR_treat_1
sample=sample_name
# path to sample fastq file
sample_r1_fq_file=/sample/fastq/R1.fastq.gz
sample_r2_fq_file=/sample/fastq/R2.fastq.gz
# path to the output folder of Cutadapt
cutadapt_out_dir=/cutadapt/output/folder
# bowtie2 rRNA index name: e.g., rRNA
mapTool_index_name_rep=rRNA
# path to bowtie2 rRNA index
mapTool_index_rep=/bowtie2/rRNA/index/folder/$hisat3n_rep_index_name
# bowtie2 small RNA index name: e.g., sRNA
mapTool_index_name_sRNA=sRNA
# path to bowtie2 small RNA index
mapTool_index_sRNA=/bowtie2/sRNA/index/folder/$mapTool_index_name_sRNA
# STAR genome index name: e.g., tair10
mapTool_index_name=tair10
# path to STAR genome index
mapTool_index=/STAR/genome/index/folder
# path to the output folder of STAR
mapTool_out_dir=/STAR/output/folder
# path to the output folder of pileup2var
pileup2var_out_dir=/pileup2var/output/folder
# path to the genome sequence file in the fasta format
genome_fa=/genome/sequence/genome.fa
# threads to run the bash script
ncpus=24
# minimum read coverage in a position to be considered for output
cov=1
# library type : F for forward, R for reverse. This parameter is used for pileup2var.
sample_pileup2var_strandness="R"


echo -e "\n$sample\n"

echo -e "\nadapter trimming\n"
if [[ ! -d $cutadapt_out_dir/$sample ]];then mkdir -p $cutadapt_out_dir/$sample;fi
# trim R1 5' end adapters and R2 3' end adapters
cutadapt -j $ncpus -e 0.1 -n 1 -O 1 -q 20 -m 29 --nextseq-trim=20  -g CTCTTCCGATCT -A AGATCGGAAGAGCGTC -o $cutadapt_out_dir/$sample/$sample.R1.a5trim.fastq.gz -p $cutadapt_out_dir/$sample/$sample.R2.a3trim.fastq.gz $sample_r1_fq_file $sample_r2_fq_file > $cutadapt_out_dir/$sample/$sample.cutadapt.round1.log
# trim R1 3' end adapters and R2 5' end adapters
cutadapt -j $ncpus -e 0.1 -n 1 -O 12 -q 20 -m 18:29 -a NNNNNNNNNNNAGATCGGAAGAGCACA -G CTCTTCCGATCT  -o $cutadapt_out_dir/$sample/$sample.R1.a3a5trim.fastq.gz -p $cutadapt_out_dir/$sample/$sample.R2.a3a5trim.fastq.gz $cutadapt_out_dir/$sample/$sample.R1.a5trim.fastq.gz $cutadapt_out_dir/$sample/$sample.R2.a3trim.fastq.gz > $cutadapt_out_dir/$sample/$sample.cutadapt.round2.log
# (optional) output clean sequences from the transcript strand (i.e., R2)
time zcat $cutadapt_out_dir/$sample/$sample.R2.a3a5trim.fastq.gz | fastx_collapser -v -Q33 -i - | fastx_trimmer -Q33 -f 12 -i - -o $cutadapt_out_dir/$sample/$sample.R2.clean.fa
wait
# extract R2 5' end barcodes (11-mer)
time umi_tools extract --random-seed=123 --extract-method=regex --bc-pattern=".*" --bc-pattern2="^(?P<umi_1>.{11}).*" -I $cutadapt_out_dir/$sample/$sample.R1.a3a5trim.fastq.gz --read2-in=$cutadapt_out_dir/$sample/$sample.R2.a3a5trim.fastq.gz \
	--stdout=$cutadapt_out_dir/$sample/$sample.R1.a3a5trim.umi.fastq.gz --read2-out=$cutadapt_out_dir/$sample/$sample.R2.a3a5trim.umi.fastq.gz -L $cutadapt_out_dir/$sample/$sample.a3a5trim.umi.log
wait
# (optional) check the distribution of R1 and R2 read length after adapter and barcode trimming
echo -e "\nDistribution of R1 length\n"
find $cutadapt_out_dir/$sample/$sample.R1.a3a5trim.umi.fastq.gz -printf "zcat %p | awk '{if(NR%%4==2) print length(\$1)}' | textHistogram -minVal=0 -maxBinCount=200 stdin \n" | sh
echo -e "\nDistribution of R2 length\n"
find $cutadapt_out_dir/$sample/$sample.R2.a3a5trim.umi.fastq.gz -printf "zcat %p | awk '{if(NR%%4==2) print length(\$1)}' | textHistogram -minVal=0 -maxBinCount=200 stdin \n" | sh

rm $cutadapt_out_dir/$sample/$sample.R1.a5trim.fastq.gz $cutadapt_out_dir/$sample/$sample.R2.a3trim.fastq.gz
rm $cutadapt_out_dir/$sample/$sample.R[12].a3a5trim.fastq.gz

if [[ -d $mapTool_out_dir/$sample ]];then rm -rf $mapTool_out_dir/$sample;fi
mkdir -p $mapTool_out_dir/$sample

echo -e "\nbowtie2 align -- rRNA mapping\n"
time bowtie2 --time --mm -p $ncpus -X 2000 --local --nofw --no-mixed --no-discordant --no-unal -x $mapTool_index_rep -1 $cutadapt_out_dir/$sample/$sample.R1.a3a5trim.umi.fastq.gz -2 $cutadapt_out_dir/$sample/$sample.R2.a3a5trim.umi.fastq.gz --un-conc-gz $mapTool_out_dir/$sample/$sample.R%.a3a5trim.umi.rmrRNA.fastq.gz | \
	samtools view -@ $ncpus -Shub /dev/stdin | samtools sort -T $mapTool_out_dir/$sample -@ $ncpus -o $mapTool_out_dir/$sample/$sample.$mapTool_index_name_rep.align.sorted.bam -

echo -e "\nbowtie2 align -- sRNA mapping\n"
time bowtie2 --time --mm -p $ncpus -X 2000 --local --nofw --no-mixed --no-discordant --no-unal -x $mapTool_index_sRNA -1 $mapTool_out_dir/$sample/$sample.R1.a3a5trim.umi.rmrRNA.fastq.gz -2 $mapTool_out_dir/$sample/$sample.R2.a3a5trim.umi.rmrRNA.fastq.gz --un-conc-gz $mapTool_out_dir/$sample/$sample.R%.a3a5trim.umi.rmrRNA.rmsRNA.fastq.gz | \
	samtools view -@ $ncpus -Shub /dev/stdin | samtools sort -T $mapTool_out_dir/$sample -@ $ncpus -o $mapTool_out_dir/$sample/$sample.$mapTool_index_name_sRNA.align.sorted.bam -

rm $mapTool_out_dir/$sample/$sample.R[12].a3a5trim.umi.rmrRNA.fastq.gz

echo -e "\nstar align -- genome mapping\n"
time STAR --genomeDir $mapTool_index --readFilesIn $mapTool_out_dir/$sample/$sample.R1.a3a5trim.umi.rmrRNA.rmsRNA.fastq.gz $mapTool_out_dir/$sample/$sample.R2.a3a5trim.umi.rmrRNA.rmsRNA.fastq.gz --readFilesCommand zcat \
    --outFileNamePrefix $mapTool_out_dir/$sample/ \
    --runThreadN $ncpus --genomeLoad NoSharedMemory     \
    --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1    \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.06              \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000         \
    --outSAMunmapped None --outFilterType BySJout --outSAMattributes All \
    --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1 \
    --limitBAMsortRAM 100000000000 --outWigType bedGraph --outWigStrand Stranded

rm $mapTool_out_dir/$sample/$sample.R[12].a3a5trim.umi.rmrRNA.rmsRNA.fastq.gz

echo -e "\nsamtools index -- raw sorted bam\n"
time samtools index $mapTool_out_dir/$sample/Aligned.sortedByCoord.out.bam

echo -e "\nsamtools flagstat -- raw sorted bam\n"
time samtools flagstat $mapTool_out_dir/$sample/Aligned.sortedByCoord.out.bam > $mapTool_out_dir/$sample/Aligned.sortedByCoord.out.flagstat.qc

echo -e "\nRemove unmapped, mate unmapped, not primary alignment, reads failing platform\n"
time samtools view -@ $ncpus -F 1804 -f 2 -Shub $mapTool_out_dir/$sample/Aligned.sortedByCoord.out.bam | samtools sort -T $mapTool_out_dir/$sample -@ $ncpus -o $mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.bam -

echo -e "\nsamtools index -- refined sorted bam\n"
samtools index $mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.bam

echo -e "\nsamtools flagstat -- refined sorted bam\n"
time samtools flagstat $mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.bam > $mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.flagstat.qc

echo -e "\nremove duplicates\n"
time umi_tools dedup --random-seed=123 --method=unique --spliced-is-unique --paired -I $mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.bam -S $mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.dedup.bam \
    --output-stats=$mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.dedup -L $mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.dedup.log

echo -e "\nsamtools index -- final sorted bam\n"
time samtools index $mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.dedup.bam

echo -e "\nsamtools flagstat -- final sorted bam\n"
time samtools flagstat $mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.dedup.bam > $mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.dedup.flagstat.qc

echo -e "\nPBC file output\n"
samtools sort -n -T $mapTool_out_dir/$sample -@ $ncpus $mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.bam | bamToBed -bedpe -i - | awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$4,$6,$9,$10}' | \
	sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}' > \
	$mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.pbc.qc

echo -e "\npileup2var\n"
if [[ ! -d $pileup2var_out_dir/$sample ]];then mkdir -m 701 -p $pileup2var_out_dir/$sample;fi
time pileup2var -t $ncpus -f 1804 -a A -c $cov -s $sample_pileup2var_strandness -g $genome_fa -b $mapTool_out_dir/$sample/$sample.$mapTool_index_name.align.sorted.dedup.bam -o $pileup2var_out_dir/$sample/$sample.$mapTool_index_name.pileup2var.txt

wait


echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
