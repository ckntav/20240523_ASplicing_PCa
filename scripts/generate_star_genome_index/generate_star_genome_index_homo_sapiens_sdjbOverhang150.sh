#!/bin/sh
#SBATCH --time=24:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=64G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL

module load star/2.7.11a

/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/star/2.7.11a/bin/STAR \
	--runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/input/star_genome_index/homo_sapiens_sdjbOverhang150 \
	--genomeFastaFiles /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/input/ensembl_fasta/homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
	--sjdbGTFfile /home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa/input/ensembl_gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf \
	--sjdbOverhang 150