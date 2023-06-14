############################################
## load libraries
library(readxl)
library(tidyverse)
library(jsonlite) # needed for api calls in HGNC functions
library(config) # needed for config loading
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/00_InhousePanels/"

## read configs
config_vars <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = "default")
config_vars_path <- config::get(file = Sys.getenv("CONFIG_FILE"),
    config = project_topic)

## set working directory
setwd(paste0(config_vars_path$projectsdir, project_name, script_path))

## set global options
options(scipen = 999)
############################################


############################################
# load global functions
# hgnc functions
source("../functions/hgnc-functions.R", local = TRUE)
############################################


############################################
# load the list of established local custom cancer panels
# TODO: read all files and combine into one table
leipzig_cancer_v6_list <- read_excel("data/Leipzig-Cancer_v6_gene-list.xlsx")
dresden_cancer_v1_list <- read_excel("data/Dresden-Cancer_v1_gene-list.xlsx")
############################################


############################################
## compute standardized gene list per panel

## 01) Leipzig cancer panel gene list
leipzig_cancer_v6_list_genes <- leipzig_cancer_v6_list %>%
  unique() %>%
  mutate(Panel = "leipzig_cancer_v6") %>%
  mutate(Source = "Leipzig-Cancer_v6_gene-list.xlsx") %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = genes, source = Source, panel = Panel)

## 02) Dresden cancer panel gene list
dresden_cancer_v1_list_genes <- dresden_cancer_v1_list %>%
  unique() %>%
  mutate(Panel = "dresden_cancer_v1_list") %>%
  mutate(Source = "Dresden-Cancer_v1_gene-list.xlsx") %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = genes, source = Source, panel = Panel)

############################################


############################################
## bind all tables and summarize
## compute count of panels a gene is reported in
## compute inhouse panel evidence as genes reported in at least 1 panel source
all_inhouse_panels_genes <- bind_rows(leipzig_cancer_v6_list_genes,
  dresden_cancer_v1_list_genes) %>%
  group_by(approved_symbol) %>%
  summarise(panel_inhouse = paste(unique(panel), collapse = "; "),
    panel_inhouse_count = n(),
    gene_name_reported = paste(unique(gene_name_reported), collapse = " | ")) %>%
  ungroup() %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(approved_symbol)) %>%
  mutate(at_least_one_panels = (panel_inhouse_count > 0))

all_inhouse_panels_genes_format <- all_inhouse_panels_genes %>%
  select(approved_symbol, hgnc_id, gene_name_reported, source = panel_inhouse, source_count = panel_inhouse_count, source_evidence = at_least_one_panels)

############################################


############################################
## save results
# TODO: gzip csv result files
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(all_inhouse_panels_genes_format,
  file = paste0("results/00_InhousePanels_genes.",
    creation_date,
    ".csv"),
  na = "NULL")
############################################
