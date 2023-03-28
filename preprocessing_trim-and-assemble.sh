#!/bin/bash
#$ -cwd -V
#$ -N assemble -j y
#$ -q short.qc
#$ -pe shmem 2

echo "****************************************************"
echo "SGE job ID: "$JOBID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "****************************************************"


# Variables
threads=2
spades=/well/shaw/users/amu125/programs/SPAdes-3.15.3-Linux/bin/spades.py
trimmomatic=/well/shaw/users/amu125/programs/Trimmomatic-0.39/trimmomatic-0.39.jar
trimmomaticAdapters=/well/shaw/users/amu125/programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa
isolate=$(sed -n "$SGE_TASK_ID"p /well/shaw/users/amu125/projects/shen/accs.txt) 
path="/well/shaw/users/amu125/projects/shen/data"

reads1=$path/"$isolate"_1.fastq.gz
reads2=$path/"$isolate"_2.fastq.gz
trimmed1Paired=$path/trimmed_"$isolate"_1.fastq.gz
trimmed1Unpaired=$path/trimmed_"$isolate"_1_unpaired.fastq.gz
trimmed2Paired=$path/trimmed_"$isolate"_2.fastq.gz
trimmed2Unpaired=$path/trimmed_"$isolate"_2_unpaired.fastq.gz

# Trim
echo "****************************************************"
echo "Started trimming reads at: "`date`
echo "****************************************************"
# 'Default parameters' as taken from the trimmomatic v0.39 manual
# TruSeq3 adapters are correct choice for HiSe
java -jar $trimmomatic PE $reads1 $reads2 $trimmed1Paired $trimmed1Unpaired $trimmed2Paired $trimmed2Unpaired ILLUMINACLIP:$trimmomaticAdapters:2:30:10 LEADING:3  TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads $threads
echo "****************************************************"
echo "Finished trimming reads at: "`date`
echo "****************************************************"

# Assembly
echo "****************************************************"
echo "Started assembling reads at: "`date`
echo "****************************************************"
python $spades -1 $trimmed1Paired -2 $trimmed2Paired --isolate -o assemblies/"$isolate" -k 21,33,55,77 -t $threads -m 1000
echo "****************************************************"
echo "Finished assembling reads at: "`date`
echo "****************************************************"

echo "****************************************************"
echo "Removing downloaded reads at: "`date`
rm $path/"$isolate"*
rm $path/trimmed*"$isolate"*
echo "****************************************************"
echo "Finished at: "`date`
echo "****************************************************"


