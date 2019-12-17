#!/bin/bash 

mkdir -p fastq/sub/

for f in fastq/*.fastq; do
	fname=$(basename $f .fastq)

	# count the number of reads
	total_lines=$(wc -l $f |cut -d" " -f 3) 
	total_reads=$((total_lines / 4))

	# array of read percentages to extract

    arr=("10" "20" "30" "40" "50" "60" "70" "80" "90" "100")
 
    # for loop that iterates over each element in arr
    for percentage in "${arr[@]}"; do
    	number_of_reads=$((total_reads * percentage / 100))
        seqtk sample $f $number_of_reads > fastq/sub/"${fname}_${percentage}_percent".fastq
    done
done 