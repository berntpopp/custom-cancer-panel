############################################
## load libraries
library(readxl)
library(tidyverse)
library(jsonlite) # needed for api calls in HGNC functions
library(config) # needed for config loading
library("R.utils")  ## gzip downloaded and result files
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/B_ManualCuration/"

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
# load the list of manually added genes
manual_gene_list <- read_excel("data/manual_gene-list.xlsx")
############################################


############################################
## compute standardized gene list and summarize
manual_gene_list_genes <- manual_gene_list %>%
  unique() %>%
  mutate(Panel = "manual") %>%
  mutate(Source = "manual_gene-list.xlsx") %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  mutate(source_count = 1) %>%
  mutate(source_evidence = TRUE) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = genes, source = Source, source_count, source_evidence)
############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(manual_gene_list_genes,
  file = paste0("results/B_ManualCuration_genes.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/B_ManualCuration_genes.", creation_date, ".csv"),
  overwrite = TRUE)
############################################
