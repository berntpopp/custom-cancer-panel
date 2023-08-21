############################################
## load libraries
library(readxl)
library(tidyverse)
library(jsonlite) # needed for api calls in HGNC functions
library(config) # needed for config loading
library("R.utils")  ## gzip downloaded and result files
library(biomaRt)    ## needed to get gene coordinates
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/V03_Ethnicity/"

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
source("../functions/ensembl-functions.R", local = TRUE)
############################################


############################################
# load the list of ethnicity snps
ethnicity_snps <- read_excel("data/ethnicity_snps.xlsx")
############################################


############################################
## generate standardized snp list
# add source column and merge the sources into one column separated by ";"
# add hg38 coordinates using snp_position_from_rs function
# then unnest the hg38_coordinates column
ethnicity_snps_panel <- ethnicity_snps %>%
  group_by(rs_id) %>%
  summarise(source = paste0(group, collapse = "; ")) %>%
  dplyr::select(snp = rs_id, source) %>%
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

write_csv(ethnicity_snps_panel,
  file = paste0("results/V03_Ethnicity.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/V03_Ethnicity.", creation_date, ".csv"),
  overwrite = TRUE)
############################################
