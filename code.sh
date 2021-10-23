#!/bin/bash

#####################################################################################################
#This Programme is to analyze the gene-expression level of targeted organism under different condition.
#In this case, it is Trypanosonma congolense. 

#Analysis Tools needed to be pre-installed: bowtie2, samtools, bedtools
#Data needed to be provided: samples RNAseq results, genome fasta, genome information bedfile
#Files generated for user:
#1.quality check result: qualitycheck_result.csv
#2.expression level comparation result: expression.csv 
#####################################################################################################

######################## PREPARATION: Download all the files needed and get current working directory
cp /localdisk/data/BPSM/AY21/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz .
cp /localdisk/data/BPSM/AY21/TriTrypDB-46_TcongolenseIL3000_2019.bed .
cp -r /localdisk/data/BPSM/AY21/fastq .
echo -e "Downloaded files to the current directory\n"
PATH=$(pwd)
unset IFS
IFS=$'\t'

######################## QUALITY CHECK PART1: Perform qualitycheck using fastqc
###### Current Directory: $PATH
echo "Performing quality check......"
mkdir qualitycheck.output
for file in ./fastq/*.fq.gz
do
	fastqc -q -t 16 -o qualitycheck.output $file
done
echo -e "Quality check finished\n"

####################### QUALITY CHECK PART2: Assess the numbers and quality of the raw sequence data
####################### Generate qualitycheck_result.csv which contains modules results of all the samples

###### Unzip the zip files produced from fastqc
for zipfile in ./qualitycheck.output/*.zip
do
	unzip -q $zipfile
	rm -f $zipfile
done
echo -e "Quality check result files unzipped\n" 

echo -e "Organizing quality check results......"
###### Save headers in file1.csv
echo -e "Paired_end_reads_ID
Basic_Statistics
Sequences_flagged_as_poor_quality
%GC
Per_base_sequence_quality
Per_sequence_quality_scores
Per_base_sequence_content
Per_sequence_GC_content
Per_base_N_content
Sequence_Length_Distribution
Sequence_Duplication_Levels
Overrepresented_Sequences
Adapter_Content"|paste -s >>file1.csv

###### Save sample IDs in file2,csv, and save quality check results in file3.csv
for dir in ./qualitycheck.output/*fastqc
do
	grep "Filename" $dir/fastqc_data.txt|awk '{FS='\t';{print $2;}}' | paste -s >> file2.csv	
	grep ">>\|GC\|Sequences flagged as poor quality" $dir/fastqc_data.txt|grep -v "END"|
	awk '{FS="\t"; if($2!="Statistics"){print $2;};if($2=="Statistics"){print $3;}}'| paste -s >>file3.csv
done

######################## QUALITY CHECK PART3:Generate the final file: qualitycheck_result.csv
paste file2.csv file3.csv > file4.csv
cat file1.csv file4.csv > qualitycheck_result.csv
rm -f file*
echo -e "Quality check results are saved in qualitycheck_result.csv. Please check.\n"

####################### ALIGNMENT: Use bowtie2 to align sample sequences with genome sequences and convert them into sam files
echo "Building genome index......"
bowtie2-build -q TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz ./fastq/genome 
echo -e "Genome index is successfully built\n"

###### Use a LONG FOR-LOOP to go through all the data
echo "Performing alignment and generating counts data......"
cd $PATH/fastq
###### FOR-LOOP BEGINS
for file in $(ls . | grep "_1.fq.gz$")
do
	name=${file:0:13}
	bowtie2 --local -p 16 -x genome -1 $name"_1.fq.gz" -2 $name"_2.fq.gz" -S $name.sam
	samtools view -bS $name.sam > $name.bam
	samtools sort $name.ban -o $name.sorted.bam
	samtools index $name.sorted.bam
####################### GENERATE COUNTS DATA: Count how many overlaps for each gene in each sample
##### Convert bam files to bed files  
	bedtools bamtobed -i $name.sorted.bam > $name.bed
##### Perform bedtools, use uniq to get counts data and save it in $name.count
	bedtools intersect -wa -wb -a $name.bed -b $PATH/TriTrypDB-46_TcongolenseIL3000_2019.bed |
	awk '{FS="\t";OFS="\t";{print $10,$11;}}' |sort| uniq -c| sed -E 's/^ *//; s/ /\t/'> $name.count
###### Extract Gene ID and Gene information from genome index bed file	
	sort -k4,4 $PATH/TriTrypDB-46_TcongolenseIL3000_2019.bed |cut -f 4,5 |
	awk '{$1=0 "\t" $1"\t"$2 }1' > index.count
###### For those genes that didn't express, append them in the $name.count and fill the blank with 0
	awk 'NR==FNR { h[$2] = $1; next }{ print $2,$1,h[$2]}' $name.count index.count |
	cut -d " " -f 3 > $name.total
	sed -i -e 's/^$/0/' $name.total
done
###### FOR-LOOP ENDS 
echo -e "Counts data generated\n"

####################### ANALYZE COUNTS DATA
echo "Analyzing counts data......"
###### Sort 100k.fqfiles by sample and time, and save it in information.txt
sort -k2,2 -k4,4 100k.fqfiles | grep -v "ID" > information.txt

###### For lines with the sample and time, save the ID in $treatment"_"$sample"_"$time.name
mkdir counts
while read ID sample replicate time treatment end1 end2
do
	echo $treatment"_"$sample"_"$time >> name
  	echo "100k."$ID".total" >> ,/counts/$treatment"_"$sample"_"$time.name
done < information.txt

###### For the three IDs in .name, calculate the mean of their counts data and save it under .mean
for i in ./counts/*.name
do
  	filename=$(cat $i | tr "\n" " ")
        paste -d+ $filename | bc | awk '{ if($1!=0) ave=$1/3;else ave=0; print ave }' > ./.counts/$i.mean
done

###### Cancatenate *.mean together, and put column names and row names in the file
paste ./counts/*.mean > ./counts/expression_level.csv
cut -f 2,3 ./fastq/index.count >./counts/rowname
paste ./counts/rowname ./counts/expression_level.csv > file1
echo -e "GeneID\nGene_Information" > file2
cat file2 name | tr "\n" "\t"> file3
cat file3 file1 > $PATH/expression.csv
rm -f file*
echo "Counts data saved in expression.csv. Please check."
