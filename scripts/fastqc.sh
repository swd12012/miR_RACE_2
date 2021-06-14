#!/bin/bash
#SBATCH --job-name=fastqc      ## Name of the job.
#SBATCH -p free          ## partition/queue name
#SBATCH --cpus-per-task=2    ## number of cores the job needs, can the programs you run make used of multiple cores?

module load fastqc

dir1='/dfs6/pub/swdu/microRNA/20210610_miR_RACE/data/raw'

#file = $dir1/filenames.txt
#prefix=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1 | cut -f2`

fastqc -o $dir1/fastqc -t 2 *.fq.gz
