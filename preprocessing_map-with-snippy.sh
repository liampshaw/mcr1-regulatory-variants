#!/bin/bash
#$ -cwd -V
#$ -N assemble -j y
#$ -q short.qc
#$ -pe shmem 2
#$ -t 1-220

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
isolate=$(sed -n "$SGE_TASK_ID"p /well/shaw/users/amu125/projects/incX4/IncX4-accs.txt)
path="/well/shaw/users/amu125/projects/incX4/raw-reads"

reads1=$path/"$isolate"_1.fastq.gz
reads2=$path/"$isolate"_2.fastq.gz

# Run snippy
# CONDA command
module load Anaconda3/2020.11
eval "$(conda shell.bash hook)"
# COMMANDS GO HERE
conda activate snippy
snippy=/well/shaw/users/amu125/programs/snippy/bin/snippy
$snippy --R1 $reads1 --R2 $reads2 --reference KU761327.1.fasta --mincov 10 --minfrac 0.9 --outdir mapping-KU761327/$isolate --cpus $threads --cleanup --force 
# N.B. cleanup removes bams



echo "Removing downloaded reads at: "`date`
# rm $path/"$isolate"*

echo "****************************************************"
echo "Finished at: "`date`
echo "****************************************************"
