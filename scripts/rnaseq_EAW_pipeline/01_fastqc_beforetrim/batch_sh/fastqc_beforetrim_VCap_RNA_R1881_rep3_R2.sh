#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/fastqc_beforetrim_output/VCap_RNA_R1881_rep3_R2


/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/fastqc/0.12.1/fastqc --outdir /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/fastqc_beforetrim_output/VCap_RNA_R1881_rep3_R2 --format fastq /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/raw/rnaseq_EAW/raw_fastq/CL208_VCaP_R18-3_S225_L004_R2_001.fastq.gz
