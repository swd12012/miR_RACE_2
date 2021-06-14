#!/bin/bash
#SBATCH --job-name=RACE_align      ## Name of the job.
#SBATCH -p free          ## partition/queue name
#SBATCH --array=1-5         ## number of tasks to launch, given hint below wc -l $file is helpful
#SBATCH --cpus-per-task=2    ## number of cores the job needs, can the programs you run make used of multiple cores?

module load samtools/1.10
module load hisat2/2.2.1

ref='/data/homezvol2/swdu/ref/mmu/Mus_musculus.GRCm38.dna.toplevel.fa'
file='/dfs6/pub/swdu/microRNA/20210331_RPES_211/samplecoding.txt'

#accidentally overwrote 
prefix=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1 | cut -f1`
rawdir='/dfs6/pub/swdu/microRNA/20210610_miR_RACE/data/raw/'
outdir='/dfs6/pub/swdu/microRNA/20210610_miR_RACE/data/processed/'

hisat2 -p 2 -x $ref -1 $rawdir${prefix}-READ1.fq.gz -2 $rawdir${prefix}-READ2.fq.gz -S $outdir${prefix}.sam
samtools view -bS $outdir${prefix}.sam > $outdir${prefix}.bam
samtools sort $outdir${prefix}.bam -o $outdir${prefix}.sorted.bam
samtools index $outdir${prefix}.sorted.bam
