setwd("/Users/chris/Desktop/20240523_ASplicing_PCa")

library(tidyverse)

#
psi_sigma_output <- read_tsv("output/rna-pipeline_EAW-GRCh38_PE/alternative_splicing_PSIsigma/ASplicing_VCap_R1881_vs_EtOH/VCap_R1881_vs_CTRL_r10_ir3.sorted.txt") %>% 
  set_names("event_region", "gene_symbol", "target_exon",
            "event_type", "N", "T",
            "exon_type", "reference_transcript", "dPSI",
            "T_test_p_value", "FDR_BH", "N_Values",
            "T_Values", "Database_ID")

#
head(psi_sigma_output)
colnames(psi_sigma_output)

psi_sigma_output %>% pull(`event_type`) %>% unique

KLK3 <- psi_sigma_output %>%
  dplyr::filter(gene_symbol == "KLK3")
KLK3

SNX22 <- psi_sigma_output %>%
  dplyr::filter(gene_symbol == "SNX22")
SNX22

psi_sigma_output %>%
  dplyr::filter(T_test_p_value <= 0.05) %>% 
  group_by(gene_symbol) %>% 
  tally %>% 
  arrange(desc(n))

psi_sigma_output %>%
  dplyr::filter(T_test_p_value <= 0.05) %>% 
  group_by(gene_symbol) %>% 
  tally %>% 
  arrange(desc(n)) %>% 
  dplyr::filter(n >= 4) %>% 
  print(n = 100)

library(Gviz)
library(ggplot2)

# Split the genomic coordinates into components
library(stringr)

# Example: Split and parse the first entry
region <- str_split(psi_sigma_output$event_region[1], ":")[[1]]
chr <- region[1]
coords <- str_split(region[2], "-")[[1]]
start <- as.numeric(coords[1])
end <- as.numeric(coords[2])

# Create a GRanges object for the event
library(GenomicRanges)

event_gr <- GRanges(seqnames = chr,
                    ranges = IRanges(start = start, end = end))

# Repeat for the target exon
exon <- str_split(psi_sigma_output$target_exon[1], ":")[[1]]
coords_exon <- str_split(exon[2], "-")[[1]]
start_exon <- as.numeric(coords_exon[1])
end_exon <- as.numeric(coords_exon[2])

exon_gr <- GRanges(seqnames = chr,
                   ranges = IRanges(start = start_exon, end = end_exon))

# Load BAM files and create a sashimi plot
library(Gviz)
library(Rsamtools)

# Path to BAM files (replace with actual paths)
bam_files <- c("/Users/chris/Desktop/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/alignments_human/Aligned.sortedByCoord.out.bam/VCap_RNA_EtOH_rep1.Aligned.sortedByCoord.out.bam")

# bam_files <- c("/Users/chris/Desktop/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/alignments_human/Aligned.sortedByCoord.out.bam/VCap_RNA_EtOH_rep1.Aligned.sortedByCoord.out.bam",
#                "/Users/chris/Desktop/20240523_ASplicing_PCa/output/rna-pipeline_EAW-GRCh38_PE/alignments_human/Aligned.sortedByCoord.out.bam/VCap_RNA_R1881_rep1.Aligned.sortedByCoord.out.bam")

# Create sashimi plot
sashimi_plot <- function(event_gr, exon_gr, bam_files) {
  # Plot the sashimi plot using Gviz
  sashimi_track <- AlignmentsTrack(range = bam_files, genome = "hg38",
                                   chromosome = start(event_gr),
                                   start = start(event_gr),
                                   end = end(event_gr))
  
  plotTracks(list(sashimi_track),
             from = start(event_gr), to = end(event_gr),
             chromosome = seqnames(event_gr),
             transcriptAnnotation = "gene",
             sashimi = TRUE,
             sashimiHeight = 0.2,
             sashimiOverlap = 0.5,
             sashimiCol = "blue",
             sashimiAlpha = 0.8)
}

# Generate the plot for the first event
sashimi_plot(event_gr, exon_gr, bam_files)

# ### Venn literature
library(vennRanges)
a <- generate_random_genomic_range()
b <- generate_random_genomic_range()
head(a)
head(b)

draw_2_venn(a, b)
draw_2_venn(a, b, "A", "B", "gray", "purple")
draw_2_venn(a, b, "A", "B", "gray", "purple", 'raw', 1)


# example not run Venn
library(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakfiles <- getSampleFiles()
peakAnnoList <- lapply(peakfiles, annotatePeak)
names(peakAnnoList) <- names(peakfiles)
genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)
