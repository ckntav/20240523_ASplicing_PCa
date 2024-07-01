setwd("/Users/chris/Desktop/20240523_ASplicing_PCa")

library(tidyverse)
library(eulerr)
library(ComplexHeatmap)

#
condition_list <- c("R1881", "R1881+E2", "E2")
output_dir <- "output/rna-pipeline_EAW-GRCh38_PE/alternative_splicing_PSIsigma"

all_psi_sigma_output <- list()
significant_psi_sigma_output <- list()
for (condition in condition_list) {
  message("# ", condition)
  condition_dir <- paste(sep = "_", "ASplicing_VCap", condition, "vs", "EtOH")
  
  output_filepath <- file.path(output_dir, condition_dir, paste(sep = "_", "VCap", condition, "vs_CTRL_r10_ir3.sorted.txt"))
  
  psi_sigma_output <- read_tsv(output_filepath, show_col_types = FALSE, progress = FALSE) %>% 
    set_names("event_region", "gene_symbol", "target_exon",
              "event_type", "N", "T",
              "exon_type", "reference_transcript", "dPSI",
              "T_test_p_value", "FDR_BH", "N_Values",
              "T_Values", "Database_ID")
  
  message(" > number of splicing events = lines : ", nrow(psi_sigma_output))
  message(" > number of significant splicing events : ", nrow(psi_sigma_output %>% dplyr::filter(T_test_p_value <= 0.05)))
  
  all_psi_sigma_output[[condition]] <- psi_sigma_output
  significant_psi_sigma_output[[condition]] <- psi_sigma_output %>% dplyr::filter(T_test_p_value <= 0.05)
}

#
R1881_output <- all_psi_sigma_output[["R1881"]]

#
colnames(R1881_output)
head(R1881_output, 8) %>% as.data.frame

R1881_output %>% pull(`event_type`) %>% unique

for (condition in condition_list) {
  message("# ", condition)
  
  psi_sigma_output <- all_psi_sigma_output[[condition]]
  
  KLK3 <- psi_sigma_output %>% dplyr::filter(gene_symbol == "KLK3") # %>% 
    # dplyr::filter(target_exon == "19:50856240-50856310")
  # print(KLK3 %>% dplyr::filter(T_test_p_value <= 0.05))
  print(KLK3 %>% as.data.frame())
}

#
all_sign_events <- list("R1881" = significant_psi_sigma_output[["R1881"]] %>% pull(Database_ID) %>% unique,
                        "R1881+E2" = significant_psi_sigma_output[["R1881+E2"]] %>% pull(Database_ID) %>% unique,
                        "E2" = significant_psi_sigma_output[["E2"]] %>% pull(Database_ID) %>% unique)

# all_sign_events <- list("R1881" = significant_psi_sigma_output[["R1881"]] %>% pull(gene_symbol) %>% unique,
#                         "R1881+E2" = significant_psi_sigma_output[["R1881+E2"]] %>% pull(gene_symbol) %>% unique,
#                         "E2" = significant_psi_sigma_output[["E2"]] %>% pull(gene_symbol) %>% unique)

sapply(all_sign_events, length)

ltm_all_sign_events <- list_to_matrix(all_sign_events)
fit_all_sign_events <- euler(ltm_all_sign_events, shape = "ellipse")
plot(fit_all_sign_events, labels = TRUE, legend = list(side = "bottom"),
     # quantities = list(type = c("counts", "percent")),
     quantities = list(type = c("counts")),
     fills = list(fill = c("#2B70AB", "#FFB027", "#3EA742", "#CD3301", "#9370DB")),
     edges = list(alpha = 0))
