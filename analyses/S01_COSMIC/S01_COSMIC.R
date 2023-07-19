############################################
## load libraries
library(tidyverse)  ##needed for general table operations
library(jsonlite)  ##needed for HGNC requests
library(config) # needed for config loading
library(rvest)
library("R.utils")  ## gzip downloaded and result files
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/S01_COSMIC/"

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
############################################


############################################
## download all required database sources from HPO and OMIM
file_date <- strftime(as.POSIXlt(Sys.time(), "UTC", "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

# COSMIC census TSV file
cosmic_census_url <- "https://cancer.sanger.ac.uk/cosmic/census/all?home=y&name=all&tier=&export=json"
cosmic_census_filename <- paste0("data/downloads/cosmic_census.", file_date, ".tsv")
download.file(cosmic_census_url, cosmic_census_filename, mode = "wb")

# COSMIC Expert Curation of Genes
cosmic_expert_url <- "https://cancer.sanger.ac.uk/cosmic/stats/point"
cosmic_expert_filename <- paste0("data/downloads/cosmic_expert.", file_date, ".html")
download.file(cosmic_expert_url, cosmic_expert_filename, mode = "wb")
############################################


############################################
## read in files, normalize and join

# load and normalize the COSMIC census table
cosmic_census_table <- read_delim(cosmic_census_filename,
    delim = "\t",
   escape_double = FALSE,
   trim_ws = TRUE)

cosmic_census_table_normalize <- cosmic_census_table %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(`Gene Symbol`)) %>%
  select(hgnc_id,
    gene_name_reported = `Gene Symbol`,
    tier = Tier,
  somatic = Somatic,
  germline = Germline,
  tumour_types_somatic = `Tumour Types(Somatic)`,
  tumour_types_germline = `Tumour Types(Germline)`,
  molecular_genetics = `Molecular Genetics`,
  role_in_cancer = `Role in Cancer`)

# load and normalize the COSMIC expert table
cosmic_expert_html <- read_html(cosmic_expert_filename)

cosmic_expert_table <- cosmic_expert_html %>%
  html_table()

cosmic_expert_table_normalize <- cosmic_expert_table[[1]] %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(Genes))

############################################


############################################
## bind both tables and summarize
## compute a score for cosmic based on:
## gene having > 500 mutations in expert list
## gene having somatic mutations
## gene having germline mutations
## gene being in COSMIC Tier 1 or 2
cosmic_gene_list <- cosmic_census_table_normalize %>%
  left_join(cosmic_expert_table_normalize, by = c("hgnc_id")) %>%
  mutate(mut_per_sample = Mutations / Samples) %>%
  mutate(over_500_mut = case_when(Mutations > 500 ~ 1,
  TRUE ~ 0)) %>%
  mutate(germline_yes = case_when(germline == "yes" ~ 1,
  TRUE ~ 0)) %>%
  mutate(somatic_yes = case_when(somatic == "yes" ~ 1,
  TRUE ~ 0)) %>%
  mutate(tier_score = case_when(tier == 1 ~ 2,
  tier == 2 ~ 1,
  TRUE ~ 0)) %>%
  mutate(source = paste0("somatic: ", tumour_types_somatic, "; ", "germline: ", tumour_types_germline)) %>%
  mutate(source_count = over_500_mut + germline_yes + somatic_yes + tier_score) %>%
  mutate(source_evidence = (source_count > 3)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol,
  hgnc_id,
  gene_name_reported,
  source,
  source_count,
  source_evidence)
############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(cosmic_gene_list,
  file = paste0("results/S01_COSMIC_genes.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/S01_COSMIC_genes.", creation_date, ".csv"),
  overwrite = TRUE)
############################################
