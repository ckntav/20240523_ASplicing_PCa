#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/alignments_tmp/VCap_RNA_R1881+E2_rep1


mkdir -p /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/alignments_human/VCap_RNA_R1881+E2_rep1


/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/star/2.7.11a/bin/STAR --runThreadN 8 --outSAMtype BAM SortedByCoordinate --outFilterIntronMotifs RemoveNoncanonical --genomeDir /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/input/star_genome_index/homo_sapiens_sdjbOverhang150 --twopassMode Basic --readFilesCommand zcat --readFilesIn /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/raw/rnaseq_EAW/fastp_output/VCap_RNA_R1881+E2_rep1_1.fastq.gz /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/raw/rnaseq_EAW/fastp_output/VCap_RNA_R1881+E2_rep1_2.fastq.gz --outTmpDir /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/alignments_tmp/VCap_RNA_R1881+E2_rep1/VCap_RNA_R1881+E2_rep1 --outFileNamePrefix /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/alignments_human/VCap_RNA_R1881+E2_rep1/VCap_RNA_R1881+E2_rep1.


/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/samtools/1.18/bin/samtools index /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/alignments_human/VCap_RNA_R1881+E2_rep1/VCap_RNA_R1881+E2_rep1.Aligned.sortedByCoord.out.bam
