############################################
## load libraries
library(tidyverse)  ## needed for general table operations
library(readr)  ## needed to read files
library(tools)  ## needed for checksums
library("R.utils")  ## gzip downloaded and result files
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/"

## read configs
config_vars_proj <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = project_topic)

## set working directory
setwd(paste0(config_vars_proj$projectsdir, project_name, script_path))

## set global options
options(scipen = 999)
############################################


############################################
## load all analyses files and transform table

# define analyses paths
analyses_paths <- c("G00_InhousePanels/results/",
  "G01_PanelApp/results/",
  "G02_HPO/results/",
  "G03_DiagnosticPanels/results/",
  "S04_COSMIC/results/",
  "A_IncidentalFindings/results/",
  "B_ManualCuration/results/")

# find all CSV files in results folders and filter
# select only genes files
# select newest file
results_csv_table <- list.files(path = analyses_paths,
    pattern = ".csv.gz",
    full.names = TRUE) %>%
  as_tibble() %>%
  separate(value, c("analysis", "path", "file"), sep = "\\/") %>%
  mutate(file_path = paste0(analysis, "/", path, "/", file)) %>%
  separate(file, c(NA, "analysis_date", NA), sep = "\\.") %>%
  mutate(results_file_id = row_number()) %>%
  mutate(md5sum_file = md5sum(file_path)) %>%
  dplyr::select(results_file_id,
    file_path,
    analysis,
    analysis_date,
    md5sum_file) %>%
  filter(str_detect(file_path, "genes")) %>%
  group_by(analysis) %>%
    filter(analysis_date == max(analysis_date)) %>%
  ungroup() %>%
  arrange(analysis)

# load the csv files
results_genes <- results_csv_table %>%
  rowwise() %>%
  mutate(genes_list = list(read_csv(file_path,
    na = "NULL",
    col_types = cols(approved_symbol = col_character(),
      hgnc_id = col_character(), gene_name_reported = col_character(),
      source = col_character(), source_count = col_double(),
      source_evidence = col_logical())
    ))) %>%
  ungroup() %>%
  select(analysis, genes_list) %>%
  unnest(genes_list) %>%
  # following is a hack to identify cancer analyses
  mutate(cancer_analysis = str_detect(analysis, "^[0-9][0-9].+"))

# generate wide table and compute
# evidence_count = sum of lists where the source_evidence is TRUE
# list_count = sum lists where gene is found (source_evidence is TRUE or FALSE)
results_genes_wider <- results_genes %>%
  select(approved_symbol, hgnc_id, analysis, source_evidence, cancer_analysis) %>%
  group_by(approved_symbol, hgnc_id) %>%
  mutate(evidence_count = sum(source_evidence == TRUE & cancer_analysis == TRUE)) %>%
  mutate(list_count =
    sum((source_evidence == TRUE | source_evidence == FALSE) & cancer_analysis == TRUE)) %>%
  ungroup %>%
  select(-cancer_analysis) %>%
  pivot_wider(
    names_from = analysis,
    values_from = source_evidence
  ) %>%
  unique() %>%
  mutate(include = (evidence_count > 1 | A_IncidentalFindings == TRUE | B_ManualCuration == TRUE))

############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(results_genes_wider,
  file = paste0("merged/CancerPanel_MergeAnalysesSources.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("merged/CancerPanel_MergeAnalysesSources.", creation_date, ".csv"),
  overwrite = TRUE)

# TODO: add source file checksums and date as output table
############################################