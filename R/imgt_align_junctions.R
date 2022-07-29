#' Insert Alignment Gaps in IMGT JUNCTION sequences
#'
#' Adjust IMGT JUNCTION sequences for length by inserting gaps according to the
#'   IMGT unique numbering rules. Vectorised for JUNCTION sequences and their
#'   gap-inserted, length-adjusted replacements.
#'
#' IMGT JUNCTION sequences are modified by gap insertions, following the IMGT
#'   unique numbering rules for JUNCTION sequences, up to the maximum length of
#'   the longest sequence. The gap insertion character used is ".".
#'   Non-canonical sequences for which the conserved 2nd-CYS and or J-PHE/J-TRP
#'   cannot be identified can either removed or left as is.
#'
#' @param .x A character vector of IMGT JUNCTION sequences defined by the
#'   conserved 2nd-CYS at position 104 through the J-PHE or J-TRP at position
#'   118.
#' @param .rm_non_canonicals Logical.
#'   * `TRUE`: removes non-canonical JUNCTION sequences that lack a 2nd-CYS
#'     and/or the J-PHE/J-TRP and replaces them with `NA`.
#'   * `FALSE`: (the default) leaves non-canonical JUNCTION sequences as they
#'     are.
#' @returns A character vector the same length as `.x`.
#' @examples
#' imgt_align_junctions(tcr_seqs$CDR3.aa)
#'
#' imgt_align_junctions(c(tcr_seqs$CDR3.aa, "LRQQLEREGRYNEQFF"))
#'
#' imgt_align_junctions(c(tcr_seqs$CDR3.aa, "LRQQLEREGRYNEQFF"),
#'                      .rm_non_canonicals = TRUE)
#'
#' try(
#'   imgt_align_junctions("ATCATATTATATATCG")
#' )
#'
#' try(
#'   imgt_align_junctions(c(104:118))
#' )
#'
#' @author Christopher Parks
#' @references
#' https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
#'
#' Lefranc M.-P., "Unique database numbering system for immunogenetic analysis",
#'   Immunology Today, 18, 509 (1997). PMID: 9386342
#'
#' Lefranc M.-P., "The IMGT unique numbering for Immunoglobulins, T cell
#'   receptors and Ig-like domains", The Immunologist, 7, 132-136 (1999).
#'
#' Lefranc, M.-P., Pommié, C., Ruiz, M., Giudicelli, V., Foulquier, E., Truong,
#'   L., Thouvenin-Contet, V. and Lefranc, G., "IMGT unique numbering for
#'   immunoglobulin and T cell receptor variable domains and Ig superfamily
#'   V-like domains", Dev. Comp. Immunol., 27, 55-77 (2003). PMID: 12477501
#'
#' Ruiz, M. and Lefranc, M.-P. "IMGT gene identification and Colliers de Perles
#'   of human immunoglobulin with known 3D structures", Immunogenetics, 53,
#'   857-883 (2002). PMID: 11862387
#'
#' Kaas, Q. and Lefranc, M.-P. "IMGT Colliers de Perles: standardized
#'   sequence-structure representations of the IgSF and MhcSF superfamily
#'   domains", Current Bioinformatics, 2, 21-30 (2007).
#'
#' Kaas, Q., Ruiz, M. and Lefranc, M.-P. "IMGT/3Dstructure-DB and
#'   IMGT/StructuralQuery, a database and a tool for immunoglobulin, T cell
#'   receptor and MHC structural data", Nucl. Acids. Res., 32, D208-D210 (2004).
#'   PMID: 14681396
#'
#' @family IMGT
#' @export
imgt_align_junctions <- function(.x, .rm_non_canonicals = FALSE) {

  x <- .x
  rm_non_canonicals <- .rm_non_canonicals

  if (!is.character(x)) {
    rlang::abort("IMGT JUNCTION sequences must be a character vector.")
  }

  is_DNA <- any(stringr::str_detect(
    x, stringr::regex("[[^ACTG]]", ignore_case = TRUE), negate = TRUE
  ))

  if (is_DNA) {
    rlang::abort("Nucleotide sequences detected. Please enter translated, amino acid sequences")
  }

  non_canonicals <- stringr::str_which(
    x, stringr::regex("^C.*F$|^C.*W$", ignore_case = TRUE), negate = TRUE
  )

  if (rm_non_canonicals) {
    if (length(non_canonicals) > 0) {
      x[non_canonicals] <- NA_character_
      rlang::inform(
        paste(
          length(non_canonicals),
          "non-canonical CDR3 sequences were removed prior to gap insertion length adjustment."
        )
      )
    } else {
      rlang::inform("No non-canonical CDR3 sequences were detected.")
    }
  }

  if(!rm_non_canonicals) {
    if (length(non_canonicals) > 0) {
      rlang::warn(
        paste(
          length(non_canonicals), "non-canonical CDR3 sequences were detected.",
          "Review them with",
          paste0("'.x[c(", paste0(non_canonicals, collapse = ", "), ")]'"),
          "(be sure to change '.x' to your input vector object).",
          "Set '.rm_non_canonicals = TRUE' to remove them if needed."
        )
      )
    }
  }

  if (max(nchar(x), na.rm = TRUE) < 15) {
    x <- c(x, "AAAAAAAAAAAAAAA")
    res <- imgt_align_junctions_internal(x, rm_non_canonicals)
    res <- res[-length(res)]
    return(res)
  }

  res <- imgt_align_junctions_internal(x, rm_non_canonicals)
  return(res)


}

imgt_align_junctions_internal <- function(x, rm_non_canonicals) {

  lengths <- max(nchar(x) - 1, na.rm = TRUE):2
  gaps <- 2:max(nchar(x) - 1, na.rm = TRUE) - 1

  replacement <- sapply(gaps, function(gaps) {
    paste0("\\1", paste0(rep(".", gaps), collapse = ""), "\\2")
  })

  n1 <- sapply(lengths, function(lengths) {
    if (lengths %% 2 != 0) (lengths + 1) / 2 - 1 else lengths / 2 - 1
  })

  n2 <- sapply(lengths, function(lengths) {
    if (lengths %% 2 != 0) (lengths + 1) / 2 - 2 else lengths / 2 - 1
  })

  aa <- paste0(
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
  )

  pattern <- mapply(function(n1, n2) {
    paste0("(^[", aa, "]{1}[", aa, "]{", n1, "})([", aa, "]{", n2, "}[", aa, "]{1}$)")
  }, n1 = n1, n2 = n2
  )

  replacement <- stats::setNames(replacement, pattern)
  x <- stringr::str_replace_all(x, stringr::regex(replacement, ignore_case = TRUE))
  return(x)
}
