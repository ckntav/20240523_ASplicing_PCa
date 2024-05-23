#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/fastp_output/VCap_RNA_R1881+E2_rep2


mkdir -p /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/raw/rnaseq_EAW/fastp_output


/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/fastp/0.23.4/bin/fastp --in1 /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/raw/rnaseq_EAW/raw_fastq/CL208_VCaP_ER-2_S227_L004_R1_001.fastq.gz --in2 /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/raw/rnaseq_EAW/raw_fastq/CL208_VCaP_ER-2_S227_L004_R2_001.fastq.gz --detect_adapter_for_pe --overrepresentation_analysis --overrepresentation_sampling 10 --thread 8 --length_required 149 --out1 /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/raw/rnaseq_EAW/fastp_output/VCap_RNA_R1881+E2_rep2_1.fastq.gz --out2 /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/raw/rnaseq_EAW/fastp_output/VCap_RNA_R1881+E2_rep2_2.fastq.gz --html /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/fastp_output/VCap_RNA_R1881+E2_rep2/VCap_RNA_R1881+E2_rep2_fastp_report.html --json /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/fastp_output/VCap_RNA_R1881+E2_rep2/VCap_RNA_R1881+E2_rep2_fastp_report.json
