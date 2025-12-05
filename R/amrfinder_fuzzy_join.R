#' Improved AMR Fuzzy Join Function
#'
#' This function performs fuzzy matching between AMRFinder Plus and RGI gene names.
#' It handles multiple corner cases and naming variations.
#'
#' @param AMRF_Element_symbol Character. Gene name from AMRFinder Plus
#' @param RGI_Best_Hit_ARO Character vector. All unique RGI gene names to search
#' @return Character. Matched RGI gene name (lowercase) or empty string if no match
#'
#' @details
#' Authors: Cowritten by Brazeau and Claude
#' Matching strategy (in order of precedence):
#' 1. Manual mapping dictionary for special cases
#' 2. Word boundary matching for short genes (3 letters)
#' 3. Parentheses normalization (remove parens and try matching)
#' 4. Slash-to-dash conversion for bifunctional proteins
#' 5. Case-insensitive substring matching (default)

amrfuzzjoin <- function(AMRF_Element_symbol, RGI_Best_Hit_ARO) {

  # ============================================================================
  # 1. Manual Mapping Dictionary for Known Special Cases
  # ============================================================================
  # These are cases that don't follow general patterns or need explicit mapping
  manual_map <- list(
    # Parentheses to no-parens variations that need exact handling
    "erm(A)" = c("erma", "erma"),
    "erm(C)" = c("ermc", "ermc"),

    # Slash vs dash in bifunctional proteins - normalize by trying both
    "aac(6')-Ie/aph(2'')-Ia" = c("aac(6')-ie-aph(2'')-ia",
                                 "aac(6')-ie/aph(2'')-ia"),

    # bla gene family - beta-lactamase variants
    # These are tricky because names can appear as:
    #   - blaZ (gene name)
    #   - PC1 beta-lactamase (blaZ) (description with gene in parens)
    #   - PC1 beta-lactamase (description without gene)
    #   - blaPC1 (alternative naming)
    "blaZ" = c("blaz", "pc1 beta-lactamase"),
    "blaPC1" = c("blapc1", "pc1 beta-lactamase", "pc1"),
    "blaI" = c("blai", "beta-lactamase repressor"),
    "blaR1" = c("blar1", "beta-lactamase sensor")
  )

  # Check manual mapping first
  if (AMRF_Element_symbol %in% names(manual_map)) {
    patterns_to_try <- manual_map[[AMRF_Element_symbol]]

    for (pttrn in patterns_to_try) {
      # Escape regex special characters
      pttrn_escaped <- stringr::str_replace_all(
        pttrn, "([\\^$.|?*+()\\[\\]{}])", "\\\\\\1")

      cond <- stringr::str_detect(
        string = stringr::str_to_lower(RGI_Best_Hit_ARO),
        pattern = pttrn_escaped
      )

      if (any(cond)) {
        out <- stringr::str_to_lower(RGI_Best_Hit_ARO[cond])
        return(ifelse(length(out) == 0, "", out[1]))
      }
    }
    # If manual mapping didn't match, continue to other strategies
  }

  # ============================================================================
  # 2. Word Boundary Matching for Short Genes (3 letters or less)
  # ============================================================================
  # This prevents "aur" from matching "aureus" in "Staphylococcus aureus"
  if (nchar(AMRF_Element_symbol) <= 3 &&
      !stringr::str_detect(AMRF_Element_symbol, "[\\(\\)\\-]")) {

    cond <- stringr::str_detect(
      string = stringr::str_to_lower(RGI_Best_Hit_ARO),
      pattern = paste0("\\b",
                       stringr::str_to_lower(AMRF_Element_symbol),
                       "\\b")
    )

    if (any(cond)) {
      out <- stringr::str_to_lower(RGI_Best_Hit_ARO[cond])
      return(ifelse(length(out) == 0, "", out[1]))
    }
  }

  # ============================================================================
  # 3. Parentheses Normalization
  # ============================================================================
  # Try removing parentheses from AMRF gene and matching
  # Example: erm(A) -> ermA
  if (stringr::str_detect(AMRF_Element_symbol, "\\(.*\\)")) {
    # Remove parentheses: erm(A) becomes ermA
    no_parens <- stringr::str_replace_all(AMRF_Element_symbol, "[\\(\\)]", "")

    # Escape other special characters if any remain
    pttrn_escaped <- stringr::str_replace_all(
      no_parens, "([\\^$.|?*+\\[\\]{}])", "\\\\\\1")

    cond <- stringr::str_detect(
      string = stringr::str_to_lower(RGI_Best_Hit_ARO),
      pattern = stringr::str_to_lower(pttrn_escaped)
    )

    if (any(cond)) {
      out <- stringr::str_to_lower(RGI_Best_Hit_ARO[cond])
      return(ifelse(length(out) == 0, "", out[1]))
    }
  }

  # ============================================================================
  # 4. Slash to Dash Conversion for Bifunctional Proteins
  # ============================================================================
  # Try converting slashes to dashes if gene contains slash
  if (stringr::str_detect(AMRF_Element_symbol, "/")) {
    slash_to_dash <- stringr::str_replace_all(AMRF_Element_symbol, "/", "-")

    # Escape regex special characters
    pttrn_escaped <- stringr::str_replace_all(
      slash_to_dash, "([\\^$.|?*+()\\[\\]{}])", "\\\\\\1")

    cond <- stringr::str_detect(
      string = stringr::str_to_lower(RGI_Best_Hit_ARO),
      pattern = stringr::str_to_lower(pttrn_escaped)
    )

    if (any(cond)) {
      out <- stringr::str_to_lower(RGI_Best_Hit_ARO[cond])
      return(ifelse(length(out) == 0, "", out[1]))
    }
  }

  # ============================================================================
  # 5. Default: Case-Insensitive Substring Matching
  # ============================================================================
  # This handles most cases where gene name is embedded in RGI description
  # Examples: fosB matches "Staphylococcus aureus FosB"
  #          blaZ matches "PC1 beta-lactamase (blaZ)"

  pttrn <- stringr::str_to_lower(AMRF_Element_symbol)

  # Escape regex special characters
  pttrn_escaped <- stringr::str_replace_all(
    pttrn, "([\\^$.|?*+()\\[\\]{}])", "\\\\\\1")

  cond <- stringr::str_detect(
    string = stringr::str_to_lower(RGI_Best_Hit_ARO),
    pattern = pttrn_escaped
  )

  out <- stringr::str_to_lower(RGI_Best_Hit_ARO[cond])
  return(ifelse(length(out) == 0, "", out[1]))
}


#' Alternative version that returns ALL matches instead of just first match
#' Useful for debugging and validation
amrfuzzjoin_all_matches <- function(AMRF_Element_symbol, RGI_Best_Hit_ARO) {

  # Use same logic as improved version but collect all matches
  all_matches <- character(0)

  # Manual mapping
  manual_map <- list(
    "erm(A)" = c("erma", "erma"),
    "erm(C)" = c("ermc", "ermc"),
    "aac(6')-Ie/aph(2'')-Ia" = c("aac(6')-ie-aph(2'')-ia",
                                 "aac(6')-ie/aph(2'')-ia"),
    "blaZ" = c("blaz", "pc1 beta-lactamase"),
    "blaPC1" = c("blapc1", "pc1 beta-lactamase", "pc1"),
    "blaI" = c("blai", "beta-lactamase repressor"),
    "blaR1" = c("blar1", "beta-lactamase sensor")
  )

  if (AMRF_Element_symbol %in% names(manual_map)) {
    patterns_to_try <- manual_map[[AMRF_Element_symbol]]
    for (pttrn in patterns_to_try) {
      pttrn_escaped <- stringr::str_replace_all(
        pttrn, "([\\^$.|?*+()\\[\\]{}])", "\\\\\\1")
      cond <- stringr::str_detect(
        string = stringr::str_to_lower(RGI_Best_Hit_ARO),
        pattern = pttrn_escaped
      )
      if (any(cond)) {
        all_matches <- c(all_matches,
                         stringr::str_to_lower(RGI_Best_Hit_ARO[cond]))
      }
    }
  }

  # Word boundary for short genes
  if (nchar(AMRF_Element_symbol) <= 3 &&
      !stringr::str_detect(AMRF_Element_symbol, "[\\(\\)\\-]")) {
    cond <- stringr::str_detect(
      string = stringr::str_to_lower(RGI_Best_Hit_ARO),
      pattern = paste0("\\b",
                       stringr::str_to_lower(AMRF_Element_symbol),
                       "\\b")
    )
    if (any(cond)) {
      all_matches <- c(all_matches,
                       stringr::str_to_lower(RGI_Best_Hit_ARO[cond]))
    }
  }

  # Parentheses normalization
  if (stringr::str_detect(AMRF_Element_symbol, "\\(.*\\)")) {
    no_parens <- stringr::str_replace_all(AMRF_Element_symbol, "[\\(\\)]", "")
    pttrn_escaped <- stringr::str_replace_all(
      no_parens, "([\\^$.|?*+\\[\\]{}])", "\\\\\\1")
    cond <- stringr::str_detect(
      string = stringr::str_to_lower(RGI_Best_Hit_ARO),
      pattern = stringr::str_to_lower(pttrn_escaped)
    )
    if (any(cond)) {
      all_matches <- c(all_matches,
                       stringr::str_to_lower(RGI_Best_Hit_ARO[cond]))
    }
  }

  # Slash to dash
  if (stringr::str_detect(AMRF_Element_symbol, "/")) {
    slash_to_dash <- stringr::str_replace_all(AMRF_Element_symbol, "/", "-")
    pttrn_escaped <- stringr::str_replace_all(
      slash_to_dash, "([\\^$.|?*+()\\[\\]{}])", "\\\\\\1")
    cond <- stringr::str_detect(
      string = stringr::str_to_lower(RGI_Best_Hit_ARO),
      pattern = stringr::str_to_lower(pttrn_escaped)
    )
    if (any(cond)) {
      all_matches <- c(all_matches,
                       stringr::str_to_lower(RGI_Best_Hit_ARO[cond]))
    }
  }

  # Default substring matching
  pttrn <- stringr::str_to_lower(AMRF_Element_symbol)
  pttrn_escaped <- stringr::str_replace_all(
    pttrn, "([\\^$.|?*+()\\[\\]{}])", "\\\\\\1")
  cond <- stringr::str_detect(
    string = stringr::str_to_lower(RGI_Best_Hit_ARO),
    pattern = pttrn_escaped
  )
  if (any(cond)) {
    all_matches <- c(all_matches,
                     stringr::str_to_lower(RGI_Best_Hit_ARO[cond]))
  }

  # Return unique matches
  all_matches <- unique(all_matches)
  return(ifelse(length(all_matches) == 0, "",
                paste(all_matches, collapse = "; ")))
}
