#' Clean Accession Numbers
#'
#' Standardizes accession numbers from plasmidfinder and mobsuite databases
#' Handles various formatting issues and edge cases
#'
#' @param acc Character vector of accession numbers
#' @param source Character. Either "plasmidfinder" or "mobsuite"
#' @return Character vector of cleaned accession numbers
#'
#' @details
#' Cleaning steps:
#' - Removes whitespace
#' - Converts to uppercase
#' - Removes version numbers (.1, .2, etc.)
#' - For mobsuite: removes cluster prefixes and contig suffixes
#' - For plasmidfinder: removes garbage prefixes
#' - Extracts valid GenBank/RefSeq patterns

clean_accessions <- function(acc, source = "plasmidfinder") {

  if (source == "plasmidfinder") {
    # Plasmidfinder cleaning
    cleaned <- acc %>%
      str_trim() %>%
      str_to_upper() %>%
      # Remove garbage prefixes like "plasmid2)"
      str_replace("^[^A-Z]*(?=[A-Z]{2})", "") %>%
      # Remove version numbers
      str_replace("\\.\\d+$", "") %>%
      # Extract valid accession pattern if buried
      str_extract("[A-Z]{1,2}[0-9]{5,}")

  } else if (source == "mobsuite") {
    # Mobsuite cleaning (applied AFTER comma separation)
    cleaned <- acc %>%
      str_trim() %>%
      str_to_upper() %>%
      # Remove MOBsuite cluster prefix (e.g., 000261__)
      str_replace("^[0-9]+__", "") %>%
      # Remove contig suffix with leading zeros (_00002, _00022)
      str_replace("_0{2,}\\d+$", "") %>%
      # Remove GenBank version (.1, .2)
      str_replace("\\.\\d+$", "") %>%
      # Remove any remaining alphanumeric prefix (e.g., NC__)
      str_replace("^[A-Z]+[0-9]*__", "") %>%
      # Extract valid accession if still malformed
      str_extract("[A-Z]{1,2}[0-9]{5,}")

  } else {
    stop("source must be 'plasmidfinder' or 'mobsuite'")
  }

  # Filter out NAs and invalid patterns
  cleaned <- ifelse(is.na(cleaned) | cleaned == "" | cleaned == "-", NA_character_, cleaned)

  return(cleaned)
}


#' Process Mobsuite Accessions
#'
#' Handles comma-separated accession lists from mobsuite
#' Separates first, then cleans individual accessions
#'
#' @param mobsuite_df Data frame with acc column
#' @return Data frame with separated and cleaned accessions
process_mobsuite_accessions <- function(mobsuite_df) {

  mobsuite_df %>%
    # Step 1: Separate comma-delimited accessions
    tidyr::separate_longer_delim(cols = acc, delim = ",") %>%
    # Step 2: Clean each accession
    dplyr::mutate(acc = clean_accessions(acc, source = "mobsuite")) %>%
    # Step 3: Filter out invalid entries
    dplyr::filter(!is.na(acc))
}
