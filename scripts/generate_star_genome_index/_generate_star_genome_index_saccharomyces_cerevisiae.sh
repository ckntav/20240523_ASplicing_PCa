#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=32G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL

module load star/2.7.9a

mkdir -p /home/chris11/projects/def-stbil30/chris11/20230914_TTseq_pipeline_demo/input/star_genome_index/saccharomyces_cerevisiae

/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/star/2.7.9a/bin/STAR \
	--runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir /home/chris11/projects/def-stbil30/chris11/20230914_TTseq_pipeline_demo/input/star_genome_index/saccharomyces_cerevisiae \
	--genomeFastaFiles /home/chris11/projects/def-stbil30/chris11/20230914_TTseq_pipeline_demo/input/ensembl_fasta/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
	--sjdbGTFfile /home/chris11/projects/def-stbil30/chris11/20230914_TTseq_pipeline_demo/input/ensembl_gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.110.gtf \
	--sjdbOverhang 100