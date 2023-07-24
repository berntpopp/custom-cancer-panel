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
script_path <- "/analyses/S00_InhousePanels/"

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
# load the list of established local custom cancer panels
neuropatho_list <- read_excel("data/Variante3_Neuropatho_eigenes_Panel_gross.xlsx",
  col_names = FALSE) %>%
  pivot_longer(everything(), names_to = "col", values_to = "gene") %>%
  select(-col) %>%
  mutate(gene = str_squish(gene)) %>%
  filter(!is.na(gene)) %>%
  filter(!str_detect(gene, " ")) %>%
  filter(!str_detect(gene, ",")) %>%
  unique()

solid_tumor_list <- read_excel("data/Genliste_RNA_Solid_Tumor_886.xlsx",
  col_names = FALSE) %>%
  pivot_longer(everything(), names_to = "col", values_to = "gene") %>%
  select(-col) %>%
  mutate(gene = str_squish(gene)) %>%
  filter(!is.na(gene)) %>%
  filter(!str_detect(gene, " ")) %>%
  filter(!str_detect(gene, ",")) %>%
  unique()

sarcoma_list <- read_excel("data/Genliste_RNA_Sarkom_Panel_856.xlsx",
  col_names = FALSE) %>%
  pivot_longer(everything(), names_to = "col", values_to = "gene") %>%
  select(-col) %>%
  mutate(gene = str_squish(gene)) %>%
  filter(!is.na(gene)) %>%
  filter(!str_detect(gene, " ")) %>%
  filter(!str_detect(gene, ",")) %>%
  unique()
############################################


############################################
## compute standardized gene list per panel

## 01) Neuropathology Leipzig cancer panel gene list
neuropatho_list <- read_excel("data/Variante3_Neuropatho_eigenes_Panel_gross.xlsx",
  col_names = FALSE) %>%
  pivot_longer(everything(), names_to = "col", values_to = "gene") %>%
  select(-col) %>%
  mutate(gene = str_squish(gene)) %>%
  filter(!is.na(gene)) %>%
  filter(!str_detect(gene, " ")) %>%
  filter(!str_detect(gene, "\\.")) %>%
  unique() %>%
  mutate(Panel = "neuropatho_leipzig") %>%
  mutate(Source = "Variante3_Neuropatho_eigenes_Panel_gross.xlsx") %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source = Source, panel = Panel)

## 02) leipzig pathology solid tumor list
solid_tumor_list <- read_excel("data/Genliste_RNA_Solid_Tumor_886.xlsx",
  col_names = FALSE) %>%
  pivot_longer(everything(), names_to = "col", values_to = "gene") %>%
  select(-col) %>%
  mutate(gene = str_squish(gene)) %>%
  filter(!is.na(gene)) %>%
  filter(!str_detect(gene, " ")) %>%
  filter(!str_detect(gene, ",")) %>%
  unique() %>%
  mutate(Panel = "solid_tumor_leipzig") %>%
  mutate(Source = "Genliste_RNA_Solid_Tumor_886.xlsx") %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source = Source, panel = Panel)

## 03) leipzig pathology sarcoma list
sarcoma_list <- read_excel("data/Genliste_RNA_Sarkom_Panel_856.xlsx",
  col_names = FALSE) %>%
  pivot_longer(everything(), names_to = "col", values_to = "gene") %>%
  select(-col) %>%
  mutate(gene = str_squish(gene)) %>%
  filter(!is.na(gene)) %>%
  filter(!str_detect(gene, " ")) %>%
  filter(!str_detect(gene, ",")) %>%
  unique() %>%
  mutate(Panel = "sarcoma_leipzig") %>%
  mutate(Source = "Genliste_RNA_Sarkom_Panel_856.xlsx") %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(gene)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = gene, source = Source, panel = Panel)

############################################


############################################
## bind all tables and summarize
## compute count of panels a gene is reported in
## compute inhouse panel evidence as genes reported in at least 1 panel source
all_inhouse_panels_genes <- bind_rows(neuropatho_list,
  solid_tumor_list,
  sarcoma_list) %>%
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
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(all_inhouse_panels_genes_format,
  file = paste0("results/S00_InhousePanels_genes.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/S00_InhousePanels_genes.", creation_date, ".csv"),
  overwrite = TRUE)
############################################