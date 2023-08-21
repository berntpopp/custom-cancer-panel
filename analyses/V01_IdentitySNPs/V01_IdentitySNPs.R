############################################
## load libraries
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
library(magrittr) # needed for list extraction
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/V01_IdentitySNPs/"

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
## download web urls

# load the list of sources
identity_panels_list <- read_excel("data/identity_panels_list.xlsx") %>%
  filter(use == "yes")

# download html using phantomJS
identity_panels_list_phantomjs <- identity_panels_list %>%
  filter(method == "phantomjs") %>%
  rowwise() %>%
  mutate(filename_download = download_url_by_phantomjs(panel_source,
    panel_name, type)) %>%
  ungroup()

# download other files using wininet method
identity_panels_list_wininet <- identity_panels_list %>%
  filter(method == "wininet") %>%
  rowwise() %>%
  mutate(downloaded = download.file(panel_source,
    paste0("data/downloads/", panel_name, ".", type),
    mode = "wb",
    quiet = TRUE,
    method = "wininet"),
    filename_download = paste0("data/downloads/", panel_name, ".", type))

# combine both lists
identity_panels <- bind_rows(identity_panels_list_phantomjs,
  identity_panels_list_wininet)
############################################


############################################
## process downloaded files

###############
# 1) pengelly_panel
pengelly_panel <- pdf_text((identity_panels %>% filter(panel_name == "pengelly_panel"))$filename_download) %>%
  read_lines() %>%
  str_squish() %>%
  as_tibble() %>%
  dplyr::select(snp = value) %>%
  filter(snp != "") %>%
  filter(str_detect(snp, pattern = "^rs[0-9]*")) %>%
  unique() %>%
  arrange(snp) %>%
  mutate(panel_name = "pengelly_panel")

###############
# 2) eurogentest_ngs_panel
eurogentest_ngs_panel <- pdf_text((identity_panels %>% filter(panel_name == "eurogentest_ngs_panel"))$filename_download) %>%
  read_lines() %>%
  str_squish() %>%
  as_tibble() %>%
  dplyr::select(snp = value) %>%
  filter(snp != "") %>%
  filter(str_detect(snp, pattern = "chr.+rs[0-9]*")) %>%
  unique() %>%
  arrange(snp) %>%
  separate(snp, into = c("chr", "pos", "ref", "alt", "snp", "maf"), sep = " ") %>%
  dplyr::select(snp) %>%
  unique() %>%
  arrange(snp) %>%
  mutate(panel_name = "eurogentest_ngs_panel")

###############
# 3) ampliseq_illumina_panel
ampliseq_illumina_panel <- pdf_text((identity_panels %>% filter(panel_name == "ampliseq_illumina_panel"))$filename_download) %>%
  read_lines() %>%
  str_squish() %>%
  as_tibble() %>%
  dplyr::select(snp = value) %>%
  filter(snp != "") %>%
  filter(str_detect(snp, pattern = "^[0-9].+rs[0-9]*")) %>%
  unique() %>%
  arrange(snp) %>%
  separate(snp, into = c("code", "chr", "start", "end", "snp", "ref", "alt"), sep = " ") %>%
  dplyr::select(snp) %>%
  unique() %>%
  arrange(snp) %>%
  mutate(panel_name = "ampliseq_illumina_panel")

###############
# 4) idt_xgen_human_id_hybridization_panel
idt_xgen_human_id_hybridization_panel <- read_tsv((identity_panels %>% filter(panel_name == "idt_xgen_human_id_hybridization_panel"))$filename_download, 
    comment = "@", col_names = FALSE) %>%
  dplyr::select(snp = X5) %>%
  filter(str_detect(snp, pattern = "^rs[0-9]*")) %>%
  dplyr::select(snp) %>%
  unique() %>%
  arrange(snp) %>%
  mutate(panel_name = "idt_xgen_human_id_hybridization_panel")

###############
# 5) idt_xgen_sample_id_amplicon_panel
idt_xgen_sample_id_amplicon_panel <- read_html((identity_panels %>% filter(panel_name == "idt_xgen_sample_id_amplicon_panel"))$filename_download) %>%
  html_nodes(xpath = '//table[contains(@class,"table-condensed")]') %>%
  html_table() %>%
  extract2(1) %>%
  dplyr::select(snp = `SNP ID`) %>%
  unique() %>%
  arrange(snp) %>%
  mutate(panel_name = "idt_xgen_sample_id_amplicon_panel")

###############
# bind all panels
# add source column and merge the sources into one column separated by ";"
# add hg38 coordinates using snp_position_from_rs function
# then unnest the hg38_coordinates column
identity_panels_all <- bind_rows(pengelly_panel,
  eurogentest_ngs_panel,
  ampliseq_illumina_panel,
  idt_xgen_human_id_hybridization_panel,
  idt_xgen_sample_id_amplicon_panel) %>%
  group_by(snp) %>%
  summarise(source = paste0(panel_name, collapse = "; ")) %>%
  arrange(snp) %>%
  dplyr::select(snp, source) %>%
  unique() %>%
  arrange(snp) %>%
  mutate(hg38_coordinates =
    snp_position_from_rs(snp)) %>%
  unnest(hg38_coordinates)
############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(identity_panels_all,
  file = paste0("results/V01_IdentitySNPs.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/V01_IdentitySNPs.", creation_date, ".csv"),
  overwrite = TRUE)
############################################
