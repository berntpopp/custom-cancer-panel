############################################
## load libraries
library(readr)
library(readxl)
library(tidyverse)
library(rvest)
library(jsonlite)
library(curl)
library(httr)
library(config) # needed for config loading
library(webdriver) # needed for headless browsing
library("R.utils")  ## gzip downloaded and result files
library(pdftools) # needed for pdf files
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/S02_CommercialPanels/"

## read configs
config_vars_proj <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = project_topic)

## set working directory
setwd(paste0(config_vars_proj$projectsdir, project_name, script_path))

## set global options
options(scipen = 999)
############################################


############################################
# load global functions
# hgnc functions
source("../functions/hgnc-functions.R", local = TRUE)
source("../functions/phantomjs-functions.R", local = TRUE)
############################################


############################################
## download all files and web urls

# load the list of sources
commercial_panels_list <- read_excel("data/somatic_commercial_panels_list.xlsx") %>%
  filter(use == "yes")

# download using curl
commercial_panels <- commercial_panels_list %>%
  filter(type != "html") %>%
  rowwise() %>%
  mutate(filename_download = paste0("data/downloads/", panel_name, ".", type)) %>%
  mutate(downloaded = download.file(panel_source,
    filename_download,
    mode = "wb",
    quiet = TRUE,
    method = method))
############################################


############################################
## compute standardized gene list per panel

## 1) idtxgen_pan_cancer_hybridization_panel
panel <- "idtxgen_pan_cancer_hybridization_panel"
filename <- (commercial_panels %>%
  filter(panel_name == panel))$filename_download

idtxgen_pan_cancer_hybridization_panel <- read_excel(filename) %>%
  select(gene = `Gene Symbol`) %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = panel) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)

## 2) illumina_trusight_oncology
panel <- "illumina_trusight_oncology"
filename <- (commercial_panels %>%
  filter(panel_name == panel))$filename_download

illumina_trusight_oncology <- read_excel(filename, skip = 2) %>%
  select(gene = `Gene symbol`) %>%
  filter(!str_detect(gene, "Small variants found in gVCF file only")) %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = panel) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)

## 3) sureselect_cancer_all_in_one_solid_tumor_assay
panel <- "sureselect_cancer_all_in_one_solid_tumor_assay"
filename <- (commercial_panels %>%
  filter(panel_name == panel))$filename_download

sureselect_cancer_all_in_one_solid_tumor_assay <- pdf_text(filename) %>%
  read_lines()

sureselect_cancer_all_in_one_solid_tumor_assay <- sureselect_cancer_all_in_one_solid_tumor_assay[grep("ABL1", sureselect_cancer_all_in_one_solid_tumor_assay):
  grep("PTPN11", sureselect_cancer_all_in_one_solid_tumor_assay)] %>%
  str_split_fixed(" {3,}", 6) %>%
  as_tibble() %>%
  select(-V1) %>%
  pivot_longer(everything(), names_to = "col", values_to = "gene") %>%
  select(-col) %>%
  filter(gene != "") %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = panel) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)

## 4) sureselect_genomics_glasgow_panel
panel <- "sureselect_genomics_glasgow_panel"
filename <- (commercial_panels %>%
  filter(panel_name == panel))$filename_download

sureselect_genomics_glasgow_panel <- pdf_text(filename) %>%
  read_lines()

sureselect_genomics_glasgow_panel_core <- sureselect_genomics_glasgow_panel[grep("AKT1", sureselect_genomics_glasgow_panel):
  grep("PMS2", sureselect_genomics_glasgow_panel)] %>%
  str_split_fixed(" {2,}", 12) %>%
  as_tibble() %>%
  select(-V1) %>%
  pivot_longer(everything(), names_to = "col", values_to = "gene") %>%
  select(-col) %>%
  mutate(gene = str_remove_all(gene, "\\*|\\+|\\•")) %>%
  mutate(gene = str_squish(gene)) %>%
  filter(gene != "") %>%
  filter(gene != "4") %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = paste0(panel, "_core")) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)

sureselect_genomics_glasgow_panel_plus <- sureselect_genomics_glasgow_panel[grep("ABL1", sureselect_genomics_glasgow_panel):
  grep("RFX5", sureselect_genomics_glasgow_panel)] %>%
  str_split_fixed(" {2,}", 12) %>%
  as_tibble() %>%
  select(-V1) %>%
  pivot_longer(everything(), names_to = "col", values_to = "gene") %>%
  select(-col) %>%
  mutate(gene = str_remove_all(gene, "\\*|\\+|\\•")) %>%
  mutate(gene = str_replace_all(gene, "TGFBRN", "TGFBR1")) %>%
  mutate(gene = str_squish(gene)) %>%
  filter(gene != "") %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = paste0(panel, "_plus")) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)

## 5) ion_ampliseq_cancer_panel
panel <- "ion_ampliseq_cancer_panel"
filename <- (commercial_panels %>%
  filter(panel_name == panel))$filename_download

ion_ampliseq_cancer_panel <- pdf_text(filename) %>%
  read_lines()

ion_ampliseq_cancer_panel <- ion_ampliseq_cancer_panel[grep("ABL1", ion_ampliseq_cancer_panel):
  grep("TP53", ion_ampliseq_cancer_panel)] %>%
  str_split_fixed(" {2,}", 8) %>%
  as_tibble() %>%
  select(-V1) %>%
  pivot_longer(everything(), names_to = "col", values_to = "gene") %>%
  select(-col) %>%
  mutate(gene = str_squish(gene)) %>%
  filter(gene != "") %>%
  separate_rows(gene, sep = " ") %>%
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = panel) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)

## 6) illumina_comprehensive_panel_v3
panel <- "illumina_comprehensive_panel_v3"
filename <- (commercial_panels %>%
  filter(panel_name == panel))$filename_download

illumina_comprehensive_panel_v3 <- read_delim(filename,
  delim = "\t", escape_double = FALSE,
  col_names = FALSE, trim_ws = TRUE) %>%
  select(gene = X4) %>%
  mutate(gene = str_remove_all(gene, "OBRA_")) %>%
  separate(gene, sep = "_", into = c("gene", NA)) %>%
  mutate(gene = str_replace_all(gene, "SP", "SP1")) %>% 
  unique() %>%
  arrange(gene) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(source = panel) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source)

############################################


############################################
## bind all tables and summarize
## compute count of panels a gene is reported in
## compute commercial panel evidence as genes reported in at least 2 panel sources
all_commercial_panels_genes <- bind_rows(idtxgen_pan_cancer_hybridization_panel,
    sureselect_cancer_all_in_one_solid_tumor_assay,
    sureselect_genomics_glasgow_panel_core,
    sureselect_genomics_glasgow_panel_plus,
    ion_ampliseq_cancer_panel,
    illumina_trusight_oncology,
    illumina_comprehensive_panel_v3
  ) %>%
  group_by(approved_symbol) %>%
  summarise(panel_commercial = paste(unique(panel), collapse = "; "),
    hgnc_id = paste(unique(hgnc_id), collapse = "; "),
    panel_commercial_count = n(),
    gene_name_reported = paste(unique(gene_name_reported), collapse = " | ")) %>%
  ungroup() %>%
  mutate(at_least_two_panels = (panel_commercial_count > 1))

all_commercial_panels_genes_format <- all_commercial_panels_genes %>%
  select(approved_symbol, hgnc_id, gene_name_reported, source = panel_commercial, source_count = panel_commercial_count, source_evidence = at_least_two_panels)
############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(all_commercial_panels_genes_format,
  file = paste0("results/S02_CommercialPanels_genes.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/S02_CommercialPanels_genes.", creation_date, ".csv"),
  overwrite = TRUE)

write_csv(commercial_panels,
  file = paste0("results/S02_CommercialPanels_list.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/S02_CommercialPanels_list.", creation_date, ".csv"),
  overwrite = TRUE)
############################################
