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
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/A_IncidentalFindings/"

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
source("../functions/phantomjs-functions.R", local = TRUE)
############################################


############################################
## download web urls

# load the list of sources
incidental_findings_source_links <- read_delim("data/incidental_findings_source_links.txt",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE)

# download using phantomJS
incidental_findings <- incidental_findings_source_links %>%
  rowwise() %>%
  mutate(filename_download = download_url_by_phantomjs(link,
    name, type)) %>%
  ungroup()
############################################


############################################
## compute standardized gene list for incidental findings

url <- (incidental_findings %>%
  filter(name == "acmg"))$link

incidental_findings_full <- read_html(url)

incidental_findings_full_genes <- incidental_findings_full %>%
  html_nodes(xpath = '//*[@id="maincontent"]//table[1]') %>%
  html_table()

incidental_findings_full_genes_normalize <- incidental_findings_full_genes[[1]] %>%
  select(genes = `Gene via GTR`) %>%
  mutate(genes = str_remove_all(genes, " .+")) %>%
  unique() %>%
  mutate(Panel = "incidental_findings") %>%
  mutate(Source = url) %>%
  mutate(hgnc_id = hgnc_id_from_symbol_grouped(genes)) %>%
  mutate(approved_symbol = symbol_from_hgnc_id_grouped(hgnc_id)) %>%
  mutate(source_count = 1) %>%
  mutate(source_evidence = TRUE) %>%
  select(approved_symbol, hgnc_id, gene_name_reported = genes, source = Source, source_count, source_evidence)
############################################


############################################
## save results
# TODO: gzip csv result files
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(incidental_findings_full_genes_normalize,
  file = paste0("results/A_IncidentalFindings_genes.",
    creation_date,
    ".csv"),
  na = "NULL")
############################################
