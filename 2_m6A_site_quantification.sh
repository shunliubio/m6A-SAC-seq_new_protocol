#!/bin/bash


######################
# Begin work section #
######################


## This script includes the following jobs:
# 1.  (optional) Merge two replicates
# 2.  Find common sites between input and treated samples. Then generate a count table
# 3.  Perform fisher and binomial tests with spikein information
# 4.  Identify m6A sites and quantification


echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)


## main setting
# path to the output of pileup2var for input replicate 1
pileup2var_res_input_rep1=/pileup2var/output/input/rep1/pileup2var.txt
# path to the output of pileup2var for input replicate 2
pileup2var_res_input_rep2=/pileup2var/output/input/rep2/pileup2var.txt
# path to the output of mergeing input replicates. Please define the file name using the suffix ".bed"
pileup2var_res_input_merged=/pileup2var/output/input/merged/pileup2var.flt.bed
# path to the output of pileup2var for treated replicate 1
pileup2var_res_treated_rep1=/pileup2var/output/treated/rep1/pileup2var.txt
# path to the output of pileup2var for treated replicate 2
pileup2var_res_treated_rep2=/pileup2var/output/treated/rep2/pileup2var.txt
# path to the output of mergeing treated replicates. Please define the file name using the suffix ".bed"
pileup2var_res_treated_merged=/pileup2var/output/treated/merged/pileup2var.flt.bed
# path to the output of finding common sites between input and treated merged samples. Please define the file name using the suffix ".txt"
input_treated_common_site_count_table=/pileup2var/output/common/count.table.txt
# path to the spikein information file including three columns with the header "motif", "bg" and "coef"
sample_spikein=/spikein/info/sample/spikein.txt
# path to the output of performing statistical tests. Please define the file name using the suffix ".txt"
input_treated_common_site_count_table_test=/pileup2var/output/common/count.table.test.txt
# path to the output of identifying m6A site and quantification. Please define the file name using the suffix ".bed"
m6A_site_quantification=pileup2var/output/common/m6A.hits.bed
# path to the chromosome size file. The column 1 is the chromosome name and column 2 is the chromosome size
chrom_size=/chromosome/size/chrom.sizes.txt
# path to the genome sequence file in the fasta format
genome_fa=/genome/sequence/genome.fa
# minimum read coverage in a position to be considered for output
cov=5
# cutoff of input mutation rate
mut_rate_input=0.05
# cutoff of input mutation count
mut_count_input=2
# cutoff of treated mutation rate
mut_rate_treated=0.02
# cutoff of treated mutation count
mut_count_treated=2
# cutoff of fisher test pvalue
fisher_test_pvalue=0.1
# cutoff of binomial test fdr
binomial_test_fdr=0.05
# cutoff of m6A site fraction
m6A_fraction=10


# (optional) Merge two replicates. Imported files are from the output of pileup2var
awk -v c=$cov -F '\t' 'BEGIN {OFS="\t"} ARGIND==1 && FNR>1 {p=$1"_"$2"_"$3;t=$6+$7+$8+$9;s[p]=t;a[p]=$6} ARGIND==2 && FNR>1 {p=$1"_"$2"_"$3;t=$6+$7+$8+$9;if(p in a) {ts=s[p]+t;ta=a[p]+$6;delete a[p]} else {ts=t;ta=$6};if(ts>=c) {m=sprintf("%.4f",ta/ts*100);print $1,$2-1,$2,p"="ta"="ts"="m,m,$3}} END {for(i in a) {if(s[i]>=c) {split(i,b,"_");m=sprintf("%.4f",a[i]/s[i]*100);print b[1],b[2]-1,b[2],i"="a[i]"="s[i]"="m,m,b[3]}}}' $pileup2var_res_input_rep1 $pileup2var_res_input_rep2 | sort -k 1,1V -k 2,2n | slopBed -b 3 -i - -g $chrom_size | fastaFromBed -s -bedOut -fi $genome_fa -bed - | awk -F '\t' 'BEGIN {OFS="\t"} length($7)==7 && $7!~/[^ACTGactg]/ {a=toupper($7);gsub(/T/,"U",a);if(a~/[ACUG][GAU][AG]AC[ACU][ACUG]/) {m="DRACH"} else if(a~/[ACUG][GAU]GAU[ACUG][ACUG]/) {m="DGAU"} else if(a~/U[AG]UA[ACU][ACU][ACUG]/) {m="URUAHH"} else if(a~/[GAUC][GAUC][GAUC][GUC][GAUC][GAUC][GAUC]/) {m="nonA"} else {m="nonDRACH"};print $1,$2+3,$3-3,$4"="a"="m,$5,$6}' > $pileup2var_res_input_merged
awk -v c=$cov -F '\t' 'BEGIN {OFS="\t"} ARGIND==1 && FNR>1 {p=$1"_"$2"_"$3;t=$6+$7+$8+$9;s[p]=t;a[p]=$6} ARGIND==2 && FNR>1 {p=$1"_"$2"_"$3;t=$6+$7+$8+$9;if(p in a) {ts=s[p]+t;ta=a[p]+$6;delete a[p]} else {ts=t;ta=$6};if(ts>=c) {m=sprintf("%.4f",ta/ts*100);print $1,$2-1,$2,p"="ta"="ts"="m,m,$3}} END {for(i in a) {if(s[i]>=c) {split(i,b,"_");m=sprintf("%.4f",a[i]/s[i]*100);print b[1],b[2]-1,b[2],i"="a[i]"="s[i]"="m,m,b[3]}}}' $pileup2var_res_treated_rep1 $pileup2var_res_treated_rep2 | sort -k 1,1V -k 2,2n | slopBed -b 3 -i - -g $chrom_size | fastaFromBed -s -bedOut -fi $genome_fa -bed - | awk -F '\t' 'BEGIN {OFS="\t"} length($7)==7 && $7!~/[^ACTGactg]/ {a=toupper($7);gsub(/T/,"U",a);if(a~/[ACUG][GAU][AG]AC[ACU][ACUG]/) {m="DRACH"} else if(a~/[ACUG][GAU]GAU[ACUG][ACUG]/) {m="DGAU"} else if(a~/U[AG]UA[ACU][ACU][ACUG]/) {m="URUAHH"} else if(a~/[GAUC][GAUC][GAUC][GUC][GAUC][GAUC][GAUC]/) {m="nonA"} else {m="nonDRACH"};print $1,$2+3,$3-3,$4"="a"="m,$5,$6}' > $pileup2var_res_treated_merged

# Find common sites between input and treated samples. Then generate a count table
# Example of the count table:
# pos          motif    type  treated_A_count  treated_depth  control_A_count  control_depth
# Chr1_7008_-  UAGAUUC  DGAU  5                5              14               14
intersectBed -wo -s -a $pileup2var_res_treated_merged -b $pileup2var_res_input_merged | awk -F '\t' 'BEGIN {OFS="\t";print "pos","motif","type","treated_A_count","treated_depth","control_A_count","control_depth"} {split($4,a,"=");split($10,b,"=");print a[1],a[5],a[6],a[2],a[3],b[2],b[3]}' > $input_treated_common_site_count_table

# Perform fisher and binomial tests with spikein information
# Example of the output:
# pos          motif    type  treated_A_count  treated_depth  control_A_count  control_depth  pvalue_f  fdr_f  pvalue_b  fdr_b
# Chr1_7008_-  UAGAUUC  DGAU  5                5              14               14             1       1    1         1
./fisher.binomial.test.R $input_treated_common_site_count_table $sample_spikein $input_treated_common_site_count_table_test

# Identify m6A sites and quantification
# Example of the output (bed format):
# Chr1  31192  31193  Chr1_31193_-=AAGAGUU=nonDRACH  52.2167  -
awk -v c=$cov -v imr=$mut_rate_input -v imc=$mut_count_input -v tmr=$mut_rate_treated -v tmc=$mut_count_treated -v fp=$fisher_test_pvalue -v bf=$binomial_test_fdr -v mf=$m6A_fraction -F '\t' 'BEGIN {OFS="\t"} ARGIND==1 && FNR>1 {if($3==0) {s[$1]=1000} else {s[$1]=1/$3}} ARGIND==2 && FNR>1 && $5>=c && $7>=c && (($7-$6)/$7<imr || $7-$6<imc) && ($5-$4>=tmc && ($5-$4)/$5>tmr) && ($8<fp || $8~/e-3[0-9][0-9]/) && ($11<bf || $11~/e-3[0-9][0-9]/) && ($5-$4)/$5*100*s[$2]>=mf && ($5-$4)/$5>=($7-$6)/$7*3 {cm=sprintf("%.4f",($7-$6)/$7*100);tm=sprintf("%.4f",($5-$4)/$5*100);f=sprintf("%.4f",($5-$4)/$5*100*s[$2]);if(f+0>100) f=100;split($1,a,"_");print a[1],a[2]-1,a[2],$1"="$2"="$3,f,a[3]}' $sample_spikein $input_treated_common_site_count_table_test > $m6A_site_quantification

wait

echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
