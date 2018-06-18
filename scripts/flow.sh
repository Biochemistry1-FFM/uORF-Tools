#!/bin/bash
#$ -N rflow
#$ -cwd 
#$ -q all.q
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -o /mnt/rflow/
#$ -e /mnt/rflow/
#$ -j y
#$ -M you@uni.edu
#$ -m a
#$ -R y
export PATH="/home/user/miniconda3/bin/:$PATH"
cd /mnt/rflow
source activate snakemake
snakemake --use-conda -s uORF-Tools/Snakefile --configfile uORF-Tools/config.yaml --directory ${PWD} -j 20 --cluster-config uORF-Tools/cluster.yaml --cluster "qsub -R y -N {cluster.jobname} -cwd -q {cluster.qname} -pe {cluster.parallelenvironment} -l {cluster.memory} -o {cluster.logoutputdir} -e {cluster.erroroutputdir} -j {cluster.joinlogs} -M you@uni.edu -m a" --latency-wait 60
