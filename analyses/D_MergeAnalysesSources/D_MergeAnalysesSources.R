############################################
## load libraries
library(tidyverse)  ## needed for general table operations
library(readr)  ## needed to read files
library(tools)  ## needed for checksums
library("R.utils")  ## gzip downloaded and result files
library(readxl)
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/D_MergeAnalysesSources/"

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
# file functions
source("../functions/file-functions.R", local = TRUE)
############################################


############################################
## load all analyses files and transform table

# define analyses paths
analyses_paths <- c("/analyses/G00_InhousePanels/results/",
  "/analyses/G01_PanelApp/results/",
  "/analyses/G02_HPO/results/",
  "/analyses/G03_DiagnosticPanels/results/",
  "/analyses/S00_InhousePanels/results/",
  "/analyses/S01_COSMIC/results/",
  "/analyses/S02_CommercialPanels/results/",
  "/analyses/A_IncidentalFindings/results/",
  "/analyses/B_ManualCuration/results/")

# find all CSV files in results folders and filter
# select only genes files
# select newest file
results_csv_table <- list.files(path = paste0(config_vars_proj$projectsdir, project_name, analyses_paths),
    pattern = ".csv.gz",
    full.names = TRUE) %>%
  as_tibble() %>%
  mutate(file_path = str_replace_all(value, paste0(config_vars_proj$projectsdir, project_name, "/analyses/"), "\\.\\.\\/")) %>%
  mutate(value = str_replace_all(value, paste0(config_vars_proj$projectsdir, project_name, "/analyses/"), "")) %>%
  separate(value, c("analysis", "path", "file"), sep = "\\/") %>%
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
  mutate(cancer_analysis = str_detect(analysis, "^[SG][0-9][0-9]_.+")) %>%
  # compute type of gene list from analysis name
  mutate(list_type =
    case_when(
      str_detect(analysis, "^[S][0-9][0-9]_.+") ~
        "somatic",
      str_detect(analysis, "^[G][0-9][0-9]_.+") ~
        "germline",
      str_detect(analysis, "^[A-z]_.+") ~
        "other",
    )
  )
############################################


############################################
# load the list of genes to be genomically covered
genomic_gene_list <- read_excel("data/genomic_gene-list.xlsx")
############################################


############################################
# germline: generate wide table and compute
# evidence_count = sum of lists where the source_evidence is TRUE
# list_count = sum lists where gene is found (source_evidence is TRUE or FALSE)
# add information if gene is included in the genomic gene list
results_genes_germline_wider <- results_genes %>%
  filter(list_type == "germline" | list_type == "other") %>%
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
  mutate(include = (evidence_count > 1 | A_IncidentalFindings == TRUE | B_ManualCuration == TRUE)) %>%
  left_join(genomic_gene_list %>%
    select(approved_symbol, include) %>%
    rename(genomic_include = include),
    by = "approved_symbol") %>%
  mutate(genomic_include = ifelse(is.na(genomic_include), FALSE, TRUE))


# somatic: generate wide table and compute
# evidence_count = sum of lists where the source_evidence is TRUE
# list_count = sum lists where gene is found (source_evidence is TRUE or FALSE)
# include information if gene is included in the germline list
results_genes_somatic_wider <- results_genes %>%
  filter(list_type == "somatic") %>%
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
  mutate(include = (evidence_count > 1)) %>%
  left_join(results_genes_germline_wider %>%
    select(approved_symbol, include) %>%
    rename(germline_include = include),
    by = "approved_symbol") %>%
  left_join(genomic_gene_list %>%
    select(approved_symbol, include) %>%
    rename(genomic_include = include),
    by = "approved_symbol") %>%
  mutate(genomic_include = ifelse(is.na(genomic_include), FALSE, TRUE)) %>%
  mutate(germline_include = ifelse(is.na(germline_include), FALSE, TRUE))

############################################


############################################
## load annotation file

# define annotation file path
annotation_path <- "../C_AnnotationHGNC/results/"

# get newest annotation file
annotation_file <- get_newest_file("non_alt_loci_set_coordinates_length", annotation_path)

# load annotation file
# remove "HGNC:" string from hgnc_id
hgnc_annotation <- read_csv(annotation_file,
    col_names = TRUE) %>%
  select(symbol, hgnc_id, hg19_genomic_size, hg38_genomic_size, hg19_cds_size_mane, hg38_cds_size_mane) %>%
  mutate(hgnc_id = str_replace_all(hgnc_id, "HGNC:", ""))

# compute panel target size germline
# use genomic size if genomic_include is TRUE
results_genes_germline_target_size <- results_genes_germline_wider %>%
  filter(include == TRUE) %>%
  left_join(hgnc_annotation,
    by = "hgnc_id") %>%
  mutate(target_size = ifelse(is.na(hg19_cds_size_mane), hg38_cds_size_mane, hg19_cds_size_mane)) %>%
  mutate(target_size = ifelse(is.na(target_size), hg19_genomic_size, target_size)) %>%
  mutate(target_size_plus_genomic = ifelse(genomic_include == TRUE, hg19_genomic_size, target_size))

# compute panel target size somatic non-germline
results_genes_somatic_target_size <- results_genes_somatic_wider %>%
  filter(include == TRUE) %>%
  filter(germline_include == FALSE) %>%
  left_join(hgnc_annotation,
    by = "hgnc_id") %>%
  mutate(target_size = ifelse(is.na(hg19_cds_size_mane), hg38_cds_size_mane, hg19_cds_size_mane)) %>%
  mutate(target_size = ifelse(is.na(target_size), hg19_genomic_size, target_size))

# sum target size cds germline
results_genes_germline_target_size_sum <- results_genes_germline_target_size %>%
  summarise(target_size = sum(target_size, na.rm = TRUE))

results_genes_germline_target_size_plus_genomic_sum <- results_genes_germline_target_size %>%
  summarise(target_size_plus_genomic = sum(target_size_plus_genomic, na.rm = TRUE))

# sum target size cds somatic non-germline
results_genes_somatic_target_size_sum <- results_genes_somatic_target_size %>%
  summarise(target_size = sum(target_size, na.rm = TRUE))
############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(results_genes_germline_wider,
  file = paste0("results/CancerPanel_Germline_MergeAnalysesSources.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/CancerPanel_Germline_MergeAnalysesSources.", creation_date, ".csv"),
  overwrite = TRUE)

write_csv(results_genes_germline_wider,
  file = paste0("results/CancerPanel_Somatic_MergeAnalysesSources.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/CancerPanel_Somatic_MergeAnalysesSources.", creation_date, ".csv"),
  overwrite = TRUE)

# TODO: add source file checksums and date as output table
############################################