setwd("/Users/chris/Desktop/20240523_ASplicing_PCa")

library(tidyverse)

##### module load star/2.7.11a
##### module load samtools/1.18
##### module load picard/3.1.0
# module load star/2.7.9a samtools/1.17 picard/2.26.3

### Align to the human genome. Sort, mark duplicates and index the genome BAM.
# cd $ALIGNDIR
# STAR --runThreadN ${THREADS} --runMode alignReads --genomeDir ${HUMANIDX} --readFilesIn ${FQ1} ${FQ2} --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outSAMunmapped None --outSAMattrRGline ID:${SAMPLE} PU:${SAMPLE} SM:${SAMPLE} LB:unknown PL:illumina --outSAMtype BAM Unsorted --outTmpDir ${TMPDIR}${SAMPLE} --outFileNamePrefix ${ALIGNDIR}${SAMPLE}.
# samtools sort --threads ${THREADS} -o ${ALIGNDIR}${SAMPLE}.sorted.bam ${ALIGNDIR}${SAMPLE}.Aligned.out.bam
# java -jar picard.jar MarkDuplicates INPUT=${ALIGNDIR}${SAMPLE}.sorted.bam OUTPUT=${ALIGNDIR}${SAMPLE}.sorted.marked.bam METRICS_FILE=${ALIGNDIR}${SAMPLE}.sorted.marked.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=LENIENT TMP_DIR=${TMPDIR}${SAMPLE}
# samtools index ${ALIGNDIR}${SAMPLE}.sorted.marked.bam
# rm ${ALIGNDIR}${SAMPLE}.sorted.bam

### Align to the yeast genome (spike-in). Sort, mark duplicates and index the genome BAM.
# cd $SPIKEDIR
# STAR --runThreadN ${THREADS} --runMode alignReads --genomeDir ${SPIKEIDX} --readFilesIn ${FQ1} ${FQ2} --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outSAMunmapped None --outSAMattrRGline ID:${SAMPLE} PU:${SAMPLE} SM:${SAMPLE} LB:unknown PL:illumina --outSAMtype BAM Unsorted --outTmpDir ${TMPDIR}${SAMPLE}.spike --outFileNamePrefix ${ALIGNDIR}${SAMPLE}.
# samtools sort --threads ${THREADS} -o ${SPIKEDIR}${SAMPLE}.sorted.bam ${SPIKEDIR}${SAMPLE}.Aligned.out.bam
# java -jar picard.jar MarkDuplicates INPUT=${SPIKEDIR}${SAMPLE}.sorted.bam OUTPUT=${SPIKEDIR}${SAMPLE}.sorted.marked.bam METRICS_FILE=${SPIKEDIR}${SAMPLE}.sorted.marked.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=LENIENT TMP_DIR=${TMPDIR}${SAMPLE}
# samtools index ${SPIKEDIR}${SAMPLE}.sorted.marked.bam 
# rm ${SPIKEDIR}${SAMPLE}.sorted.bam

#
fastq_list_filename <- "rnaseq_EAW_fastq_list.txt"
df <- read_tsv(file.path("input", "rnaseq_EAW", fastq_list_filename))
fastq_folder <- "rnaseq_EAW"
output_pipeline_dir <- "rna-pipeline_EAW-GRCh38_PE"
script_pipeline_dir <- "rnaseq_EAW_pipeline"
workdir <- "/home/chris11/projects/def-stbil30/chris11/20240523_ASplicing_PCa"

#
header_sh <- c("#!/bin/sh",
               "#SBATCH --time=12:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --ntasks-per-node=1",
               "#SBATCH --cpus-per-task=4",
               "#SBATCH --mem-per-cpu=32G",
               "#SBATCH --account=def-stbil30",
               "#SBATCH --mail-user=christophe.tav@gmail.com",
               "#SBATCH --mail-type=ALL")

# Tools
# fastqc_path <- "/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/fastqc/0.11.9/fastqc"
# STAR <- "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/star/2.7.9a/bin/STAR"
STAR <- "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/star/2.7.11a/bin/STAR"
# samtools <- "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/samtools/1.17/bin/samtools"
samtools <- "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/samtools/1.18/bin/samtools"
# picard <- "java -jar /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/picard/2.26.3/picard.jar"
picard <- "java -jar /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Core/picard/3.1.0/picard.jar"
# Threads - set to take advantage of multi-threading and speed things up.
nThreads <- 4

# Path to STAR genome indices
HUMANIDX <- file.path(workdir, "input/star_genome_index/homo_sapiens_sdjbOverhang150")
# SPIKEIDX <- "generate_star_genome_index_saccharomyces_cerevisiae.sh"

for (i in 1:nrow(df)) {
  sample_name <- df$sample_name[i]
  # message("# ", i, " | ", sample_name)
  
  #
  output_fastp_dir <- file.path(workdir, "raw", fastq_folder, "fastp_output")
  fastq_R1_filename <- paste0(sample_name, "_1.fastq.gz")
  fastq_R1_filepath <- file.path(output_fastp_dir, fastq_R1_filename)
  fastq_R2_filename <- paste0(sample_name, "_2.fastq.gz")
  fastq_R2_filepath <- file.path(output_fastp_dir, fastq_R2_filename)
  # message(" > fastq R1 : ", fastq_R1_filename)
  # message(" > fastq R2 : ", fastq_R2_filename)
  
  #
  align_tmp_dir <- file.path(workdir, "output", output_pipeline_dir, "alignments_tmp", sample_name)
  align_human_dir <- file.path(workdir, "output", output_pipeline_dir, "alignments_human", sample_name)
  # align_spike_dir <- file.path(workdir, "output", "ttseq-pipeline-demo-GRCh38_PE", "alignments_spike", sample_name)
  
  #
  call_mkdir_align_tmp <- paste("mkdir", "-p", align_tmp_dir)
  call_mkdir_align_human <- paste("mkdir", "-p", align_human_dir)
  # call_mkdir_align_spike <- paste("mkdir", "-p", align_spike_dir)
  
  ##########################################################################################
  # Align to the human genome. Sort, mark duplicates and index the genome BAM.
  ##########################################################################################
  call_star_human <- paste(STAR,
                           "--runThreadN", nThreads,
                           "--outSAMtype", "BAM", "SortedByCoordinate",
                           "--outFilterIntronMotifs", "RemoveNoncanonical",
                           "--genomeDir", HUMANIDX,
                           "--twopassMode", "Basic",
                           "--readFilesCommand", "zcat",
                           "--readFilesIn", fastq_R1_filepath, fastq_R2_filepath,
                           "--outTmpDir", file.path(align_tmp_dir, sample_name),
                           "--outFileNamePrefix", file.path(align_human_dir, paste0(sample_name, ".")))
  
  # call_star_human <- paste(STAR,
  #                          "--runThreadN", nThreads,
  #                          "--runMode", "alignReads",
  #                          "--genomeDir", HUMANIDX,
  #                          "--readFilesIn", fastq_R1_filepath, fastq_R2_filepath,
  #                          "--readFilesCommand", "zcat",
  #                          "--quantMode", "TranscriptomeSAM", "GeneCounts",
  #                          "--twopassMode", "Basic",
  #                          "--outSAMunmapped", "None",
  #                          "--outSAMattrRGline", paste0("ID:", sample_name), paste0("PU:", sample_name), paste0("SM:", sample_name), "LB:unknown", "PL:illumina",
  #                          "--outSAMtype", "BAM", "Unsorted",
  #                          # "--outTmpDir", align_tmp_dir,
  #                          "--outTmpDir", file.path(align_tmp_dir, sample_name),
  #                          "--outFileNamePrefix", file.path(align_human_dir, paste0(sample_name, ".")))
  
  # call_samtools_sort_human <- paste(samtools, "sort",
  #                                   "--threads", nThreads,
  #                                   "-o", file.path(align_human_dir, paste0(sample_name, ".sorted.bam")),
  #                                   file.path(align_human_dir, paste0(sample_name, ".Aligned.out.bam")))
  # 
  # call_markdup_human <- paste(picard, "MarkDuplicates",
  #                             paste0("INPUT=", file.path(align_human_dir, paste0(sample_name, ".sorted.bam"))),
  #                             paste0("OUTPUT=", file.path(align_human_dir, paste0(sample_name, ".sorted.marked.bam"))),
  #                             paste0("METRICS_FILE=", file.path(align_human_dir, paste0(sample_name, ".sorted.marked.metrics"))),
  #                             "REMOVE_DUPLICATES=false",
  #                             "ASSUME_SORTED=true",
  #                             "MAX_RECORDS_IN_RAM=2000000",
  #                             "VALIDATION_STRINGENCY=LENIENT",
  #                             paste0("TMP_DIR=", align_tmp_dir))
  
  call_samtools_index_human <- paste(samtools, "index",
                                     file.path(align_human_dir, paste0(sample_name, ".Aligned.sortedByCoord.out.bam")))
  
  # call_rm_human <- paste("rm", file.path(align_human_dir, paste0(sample_name, ".sorted.bam")))
  
  ##########################################################################################
  # Generate script
  ########################################################################################## 
  
  file_sh <- file.path("scripts", script_pipeline_dir, "04_align_with_STAR_human/batch_sh",
                       paste0("align_with_STAR_human_", sample_name, ".sh"))
  message("sbatch ", file_sh)
  fileConn <- file(file_sh)
  writeLines(c(header_sh, "\n",
               call_mkdir_align_tmp, "\n", call_mkdir_align_human, "\n",
               call_star_human, "\n",
               # call_samtools_sort_human, "\n",
               # call_markdup_human, "\n",
               call_samtools_index_human),
             fileConn)
  close(fileConn)
}


