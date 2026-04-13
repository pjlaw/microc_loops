#!/bin/bash
#SBATCH -J microc_pipeline
#SBATCH -t 120:00:00
#SBATCH -p master-worker
#SBATCH -c 1
#SBATCH -o /home/plaw/scripts/tmp/microc_pipeline_%j.o

module load Nextflow

# Limit NF driver 
export NXF_OPTS="-Xms8G -Xmx8G"

cd /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/nextflow/pipelines/microc_loops/

nextflow run main.nf -profile cluster --input "/data/scratch/DGE/DUDGE/MOPOPGEN/mmandelia/hic_convert/cool/*.mcool" --outdir  /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/microC/renal --resolutions 1000,5000,10000 -resume

#nextflow run main_no_raichu.nf -profile cluster --input /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/microC/renal/786-O_Mega_inter_30_cool/raichu/786-O_Mega_inter_30_cool.mcool --outdir  /data/scratch/DGE/DUDGE/MOPOPGEN/plaw/microC/renal --resolutions 1000,5000,10000 -resume

