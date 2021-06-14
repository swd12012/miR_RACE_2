#!/bin/bash
#SBATCH --job-name=count      ## Name of the job.
#SBATCH -A kpalczew_lab         		 ## account to charge 
#SBATCH -p standard           		 ## partition/queue name
#SBATCH --cpus-per-task=4  		     ## number of cores the job needs, can the programs you run make used of multiple cores?

#Run this script from the directory containing the bam files
#Amend this script in the future with an absolute path to bam files

Dir1='/dfs6/pub/swdu/microRNA/20210610_miR_RACE'

module load subread/2.0.1
gtf='/data/homezvol2/swdu/ref/mm39/mm39.ncbiRefSeq.gtf.gz'
myfile=`cat /dfs6/pub/swdu/microRNA/20210610_miR_RACE/samplecoding.txt | cut -f3`
featureCounts -p -T 4 -t exon -g gene_id -Q 30 -F GTF -a $gtf -o /dfs6/pub/swdu/microRNA/20210610_miR_RACE/data/out/counts.txt $myfile
