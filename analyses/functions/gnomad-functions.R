require(httr)
require(jsonlite)
require(tidyverse)
require(purrr)

#' Fetch Gene Data from gnomAD GraphQL API
#'
#' This function takes an Ensemble gene identifier as input, makes a POST request
#' to the gnomAD GraphQL API, and returns the gene data as a tibble.
#'
#' @param ensemble_id A character string representing the Ensemble gene identifier.
#'
#' @return A tibble with columns for each field returned by the API. This includes
#'         fields such as gene_id, gene_version, symbol, and others, as well as nested
#'         fields under gnomad_constraint.
#'
#' @examples
#' \dontrun{
#'   gene_data <- get_gene_data_from_gnomad("ENSG00000008710")
#'   print(gene_data)
#' }
#'
#' @export
get_gene_data_from_gnomad <- function(ensemble_id) {

  # API URL
  api_url <- "https://gnomad.broadinstitute.org/api/"

  # Request Body
  body <- list(
    query = paste0(
      '{
        gene(gene_id: "', ensemble_id, '", reference_genome: GRCh37) {
          reference_genome
          gene_id
          gene_version
          symbol
          gencode_symbol
          name
          canonical_transcript_id
          hgnc_id
          ncbi_id
          omim_id
          chrom
          start
          stop
          strand
          gnomad_constraint {
            exp_lof
            exp_mis
            exp_syn
            obs_lof
            obs_mis
            obs_syn
            oe_lof
            oe_lof_lower
            oe_lof_upper
            oe_mis
            oe_mis_lower
            oe_mis_upper
            oe_syn
            oe_syn_lower
            oe_syn_upper
            lof_z
            mis_z
            syn_z
            pLI
          }
        }
      }'
    )
  )

  # Send POST request
  response <- POST(api_url, body = body, encode = "json")

  # Check if request was successful
  if (response$status_code == 200) {
    # Parse the JSON response
    data <- fromJSON(content(response, "text", encoding = "UTF-8"))

    # Convert the nested list to a flat data frame
    df <- as.data.frame(data$data$gene)

    # Convert the data frame to a tibble and return
    return(as_tibble(df))
  } else {
    stop("Request failed with status code ", response$status_code)
  }
}


#' Fetch ClinVar Variants Data from gnomAD GraphQL API
#'
#' This function takes an Ensemble gene identifier as input, makes a POST request
#' to the gnomAD GraphQL API, and returns the ClinVar variants data as a tibble.
#'
#' @param ensemble_id A character string representing the Ensemble gene identifier.
#'
#' @return A tibble with columns for each field returned by the API. This includes
#'         fields such as clinical_significance, clinvar_variation_id, and others, as well as nested
#'         fields under gnomad and exome/genome.
#'
#' @examples
#' \dontrun{
#'   variants_data <- getClinVarVariants("ENSG00000008710")
#'   print(variants_data)
#' }
#'
#' @export
getClinVarVariants <- function(ensemble_id) {

  # API URL
  api_url <- "https://gnomad.broadinstitute.org/api/"

  # Request Body
  body <- list(
    query = paste0(
      '{
        gene(gene_id: "', ensemble_id,'", reference_genome: GRCh37) {
          clinvar_variants {
            clinical_significance
            clinvar_variation_id
            gold_stars
            hgvsc
            hgvsp
            in_gnomad
            major_consequence
            pos
            review_status
            transcript_id
            variant_id
          }
        }
      }'
    )
  )

  # Send POST request
  response <- POST(api_url, body = body, encode = "json")

  # Check if request was successful
  if (response$status_code == 200) {
    # Parse the JSON response
    data <- fromJSON(content(response, "text", encoding = "UTF-8"), flatten = TRUE)

    # Convert the nested list to a flat data frame
    df <- as.data.frame(data$data$gene)

    # Convert the data frame to a tibble and return
    return(as_tibble(df))
  } else {
    stop("Request failed with status code ", response$status_code)
  }
}


#' Fetch Gene Data from gnomAD GraphQL API for multiple genes
#'
#' This function takes a vector of Ensemble gene identifiers as input,
#' makes a POST request for each gene identifier
#' to the gnomAD GraphQL API, and returns the gene data as a tibble.
#'
#' @param ensemble_ids A character vector representing the Ensemble gene identifiers.
#'
#' @return A tibble with rows for each gene and columns for each field returned by the API.
#'         This includes fields such as gene_id, gene_version, symbol, and others, 
#'         as well as nested fields under gnomad_constraint.
#'
#' @examples
#' \dontrun{
#'   gene_data <- get_multiple_gene_data_from_gnomad(c("ENSG00000008710", "ENSG00000012048"))
#'   print(gene_data)
#' }
#'
#' @export
get_multiple_gene_data_from_gnomad <- function(ensemble_ids) {
  # Use purrr::map to iterate over ensemble_ids and get a list of tibbles
  list_of_tibbles <- purrr::map(ensemble_ids, get_gene_data_from_gnomad)
  
  # Use purrr::reduce to bind all tibbles into a single tibble
  gene_data <- purrr::reduce(list_of_tibbles, dplyr::bind_rows)

  # Return the combined tibble
  return(gene_data)
}


#' Fetch Variant rsids Data from gnomAD GraphQL API
#'
#' This function takes a variant in VCF format and a dataset identifier as input,
#' makes a POST request to the gnomAD GraphQL API, and returns the rsids for the variant as a tibble.
#'
#' @param vcf_variant A character string representing the variant in VCF format (e.g., "1-55516888-G-GA").
#' @param dataset A character string representing the dataset identifier (e.g., "gnomad_r2_1").
#'
#' @return A tibble with columns for each field returned by the API. This includes
#'         fields such as variant_id and rsids.
#'
#' @examples
#' \dontrun{
#'   variant_data <- get_variant_rsids_from_gnomad("1-55516888-G-GA", "gnomad_r2_1")
#'   print(variant_data)
#' }
#'
#' @export
get_variant_rsids_from_gnomad <- function(vcf_variant, dataset) {

  # API URL
  api_url <- "https://gnomad.broadinstitute.org/api/"

  # Request Body
  body <- list(
    query = paste0(
      '{
        variant(variantId: "', vcf_variant, '", dataset: ', dataset, ') {
          variant_id
          rsids
        }
      }'
    )
  )

  # Send POST request
  response <- POST(api_url, body = body, encode = "json")

  # Check if request was successful
  if (response$status_code == 200) {
    # Parse the JSON response
    data <- fromJSON(content(response, "text", encoding = "UTF-8"))

    # Convert the nested list to a flat data frame
    df <- as.data.frame(data$data$variant)

    # Convert the data frame to a tibble and return
    return(as_tibble(df))
  } else {
    stop("Request failed with status code ", response$status_code)
  }
}


#' Fetch Variant rsids Data from gnomAD GraphQL API for multiple variants
#'
#' This function takes a vector of variants in VCF format and a dataset identifier as input,
#' makes a POST request for each variant
#' to the gnomAD GraphQL API, and returns the rsids for the variants as a tibble.
#'
#' @param vcf_variants A character vector representing the variants in VCF format.
#' @param dataset A character string representing the dataset identifier (e.g., "gnomad_r2_1").
#'
#' @return A tibble with rows for each variant and columns for each field returned by the API.
#'         This includes fields such as variant_id and rsids.
#'
#' @examples
#' \dontrun{
#'   variant_data <- get_multiple_variant_rsids_from_gnomad(c("1-55516888-G-GA", "1-204502514-T-TTCTGAAACAGGG"), "gnomad_r2_1")
#'   print(variant_data)
#' }
#'
#' @export
get_multiple_variant_rsids_from_gnomad <- function(vcf_variants, dataset) {
  # Use purrr::map to iterate over vcf_variants and get a list of tibbles
  list_of_tibbles <- purrr::map(vcf_variants, function(variant) {
    get_variant_rsids_from_gnomad(variant, dataset)
  })

  # Use purrr::reduce to bind all tibbles into a single tibble
  variant_data <- purrr::reduce(list_of_tibbles, dplyr::bind_rows)

  # Return the combined tibble
  return(variant_data)
}


#' Fetch Lifted Over Variant Position from gnomAD GraphQL API
#'
#' This function takes a variant in VCF format and a reference genome as input,
#' makes a POST request to the gnomAD GraphQL API, and returns the lifted over position for the variant as a tibble.
#'
#' @param vcf_variant A character string representing the variant in VCF format (e.g., "1-55516888-G-GA").
#' @param reference_genome A character string representing the reference genome (e.g., "GRCh37").
#'
#' @return A tibble with columns for each field returned by the API. This includes
#'         fields such as variant_id and reference_genome.
#'
#' @examples
#' \dontrun{
#'   liftover_data <- get_variant_liftover_from_gnomad("1-55516888-G-GA", "GRCh37")
#'   print(liftover_data)
#' }
#'
#' @export
get_variant_liftover_from_gnomad <- function(vcf_variant, reference_genome) {

  # API URL
  api_url <- "https://gnomad.broadinstitute.org/api/"

  # Request Body
  body <- list(
    query = paste0(
      '{
        liftover(source_variant_id: "', vcf_variant, '", reference_genome: ', reference_genome, ') {
          liftover {
            variant_id
            reference_genome
          }
        }
      }'
    )
  )

  # Send POST request
  response <- POST(api_url, body = body, encode = "json")

  # Check if request was successful
  if (response$status_code == 200) {
    # Parse the JSON response
    data <- fromJSON(content(response, "text", encoding = "UTF-8"))

    # Convert the nested list to a flat data frame
    df <- as.data.frame(data$data$liftover$liftover)

    # Convert the data frame to a tibble and return
    return(as_tibble(df))
  } else {
    stop("Request failed with status code ", response$status_code)
  }
}


#' Fetch Lifted Over Variant Position from gnomAD GraphQL API for multiple variants
#'
#' This function takes a vector of variants in VCF format and a reference genome as input,
#' makes a POST request for each variant
#' to the gnomAD GraphQL API, and returns the lifted over positions for the variants as a tibble.
#'
#' @param vcf_variants A character vector representing the variants in VCF format.
#' @param reference_genome A character string representing the reference genome (e.g., "GRCh37").
#'
#' @return A tibble with rows for each variant and columns for each field returned by the API.
#'         This includes fields such as variant_id and reference_genome.
#'
#' @examples
#' \dontrun{
#'   liftover_data <- get_multiple_variant_liftover_from_gnomad(c("1-55516888-G-GA", "1-204502514-T-TTCTGAAACAGGG"), "GRCh37")
#'   print(liftover_data)
#' }
#'
#' @export
get_multiple_variant_liftover_from_gnomad <- function(vcf_variants, reference_genome) {
  # Use purrr::map to iterate over vcf_variants and get a list of tibbles
  list_of_tibbles <- purrr::map(vcf_variants, function(variant) {
    get_variant_liftover_from_gnomad(variant, reference_genome)
  })
  
  # Use purrr::reduce to bind all tibbles into a single tibble
  liftover_data <- purrr::reduce(list_of_tibbles, dplyr::bind_rows)

  # Return the combined tibble
  return(liftover_data)
}
