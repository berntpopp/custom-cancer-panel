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
library(biomaRt)    ## needed to get gene coordinates
############################################


############################################
## define relative script path
project_topic <- "lb"
project_name <- "custom-cancer-panel"
script_path <- "/analyses/V02_PRS/"

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
source("../functions/phantomjs-functions.R", local = TRUE)
source("../functions/ensembl-functions.R", local = TRUE)
source("../functions/gnomad-functions.R", local = TRUE)
############################################


############################################
## download web urls

# load the list of sources
prs_list <- read_excel("data/prs_list.xlsx") %>%
  filter(use == "yes")

# download files using wininet method
prs_list <- prs_list %>%
  filter(method == "wininet") %>%
  rowwise() %>%
  mutate(downloaded = download.file(panel_source,
    paste0("data/downloads/", panel_name, ".", type),
    mode = "wb",
    quiet = TRUE,
    method = "wininet"),
    filename_download = paste0("data/downloads/", panel_name, ".", type))

############################################


############################################
## process downloaded files
# add vcf nomenclature
# ad rs_id using snp_position_from_hg19 function
# add hg38 coordinates using snp_position_from_rs function
###############
# 1) bcac_313
bcac_313 <- read_csv((prs_list %>% filter(panel_name == "bcac_313"))$filename_download, 
    skip = 1) %>%
  mutate(prs_name = "bcac_313") %>%
  mutate(hg19_vcf = paste0(Chromosome, "-",
    Position, "-",
    Reference_Allele, "-",
    Effect_Allele)) %>%
  head(50) %>% # for testing
  rowwise() %>%
  mutate(snp = get_multiple_variant_rsids_from_gnomad(hg19_vcf, dataset = "gnomad_r2_1"))

# preliminary results
bcac_313_position <- bcac_313 %>%
  unnest(c(snp),
    names_sep = "_") %>%
  mutate(hg19 =
      snp_position_from_rs(snp_rsids, reference = "hg19"),
    hg38 =
      snp_position_from_rs(snp_rsids)) %>%
  unnest(c(hg19, hg38),
    names_sep = "_")

creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(bcac_313_position,
  file = paste0("results/bcac_313_position.",
    creation_date,
    ".csv"),
  na = "NULL")
# preliminary results

###############
# bind all prs lists
# add source column and merge the sources into one column separated by ";"
# add hg19/hg38 coordinates using snp_position_from_rs function
# then unnest the coordinates columns
prs_all <- bind_rows(bcac_313) %>%
  group_by(snp) %>%
  summarise(source = paste0(panel_name, collapse = "; ")) %>%
  arrange(snp) %>%
  dplyr::select(snp, source) %>%
  unique() %>%

############################################


############################################
## save results
creation_date <- strftime(as.POSIXlt(Sys.time(),
  "UTC",
  "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d")

write_csv(prs_all,
  file = paste0("results/V02_PRS.",
    creation_date,
    ".csv"),
  na = "NULL")

gzip(paste0("results/V02_PRS.", creation_date, ".csv"),
  overwrite = TRUE)
############################################
