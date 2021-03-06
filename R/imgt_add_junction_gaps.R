#' Insert Length Adjustment Gaps in IMGT JUNCTION sequences
#'
#' Adjust IMGT JUNCTION sequences for length by inserting gaps according to the
#'   IMGT unique numbering rules. Vectorised for JUNCTION sequences and their
#'   gap-inserted, length-adjusted replacements.
#'
#' IMGT JUNCTION sequences shorter than 15 amino acids are modified by gap
#'   insertions, following the IMGT unique numbering rules for JUNCTION
#'   sequences, to a total length of 15 positions. The gap insertion character
#'   used is ".". Non-canonical sequences for which the conserved 2nd-CYS and or
#'   J-PHE/J-TRP cannot be identified can either removed or left as is without
#'   being adjusted. It is recommended to first filter these sequences out if
#'   downstream analysis depends on comparing Junction sequences that have been
#'   length adjusted using this function (see examples). IMGT unique numbering
#'   for JUNCTION sequences accommodates sequences longer than 15 amino acids,
#'   but that is not implemented here yet. JUNCTION sequences longer than 15
#'   positions are left as is.
#'
#' @param .x A character vector of IMGT JUNCTION sequences defined by the
#'   conserved 2nd-CYS at position 104 through the J-PHE or J-TRP at position
#'   118.
#' @param .rm_non_canonicals Logical.
#'   * `TRUE`: removes non-canonical JUNCTION sequences that lack a 2nd-CYS
#'     and/or the J-PHE/J-TRP and replaces them with `NA`.
#'   * `FALSE`: (the default) leaves non-canonical JUNCTION sequences as they
#'     are.
#'
#' @returns A character vector the same length as `.x`.
#'
#' @examples
#' # Simple case of 1 sequence:
#' imgt_add_junction_gaps("CASSSGTAPFSFW")
#'
#' # Multiple sequences containing short, long, and non-canonical JUNCTIONS:
#' cdr3_seqs <- c(
#'   "CASSF",
#'   "CASSGEKLFF",
#'   "CASSKPDRGIYGYTF",
#'   "TGPLHF",
#'   "CASSQETRYDFLTIDTGGKKKNTEAFF"
#' )
#' imgt_add_junction_gaps(cdr3_seqs)
#'
#' # Remove non-canonical JUNCTIONS:
#' imgt_add_junction_gaps(cdr3_seqs, .rm_non_canonicals = TRUE)
#'
#' # Input sequences must be amino acids:
#' try(
#'   imgt_add_junction_gaps("ATCATATTATATATCG")
#' )
#'
#' # Input sequences must be a character vector:
#' try(
#'   imgt_add_junction_gaps(c(104:118))
#' )
#'
#' @author Christopher Parks
#'
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
#' @seealso {imgt_simple_headers()} {imgt_get_ref_seqs()}
#'
#' @export
#'
imgt_add_junction_gaps <- function(.x, .rm_non_canonicals = FALSE) {
  if (!is.character(.x)) {
    rlang::abort("IMGT JUNCTION sequences must be a character vector.")
  }

  .is_DNA <- any(
    stringr::str_detect(
      .x,
      stringr::regex("[[^ACTG]]", ignore_case = TRUE),
      negate = TRUE
    )
  )

  if (.is_DNA) {
    rlang::abort("Nucleotide sequences detected. Please enter translated, amino acid sequences")
  }

  .bad_cdr3s <- stringr::str_which(
    .x,
    stringr::regex("^C.*W$|^C.*F$", ignore_case = TRUE),
    negate = TRUE
  )

  if (.rm_non_canonicals) {
    if (length(.bad_cdr3s) != 0) {
      .x[.bad_cdr3s] <- NA_character_
      rlang::inform(
        paste(
          length(.bad_cdr3s),
          "non-canonical CDR3 sequences were removed prior to gap insertion length adjustment."
        )
      )
    } else {
      rlang::inform("No non-canonical CDR3 sequences were detected.")
    }
  } else {
    if (length(.bad_cdr3s) != 0) {
      rlang::warn(
        paste0(
          length(.bad_cdr3s), " ", "non-canonical Junctions were detected.\n",
          "  These have been left as is, without length adjustments.\n",
          "  set remove_non_canonicals = TRUE to remove them."
        )
      )
    }
  }

  .seqs <- .x

  .gaps <- 1:13

  .replacement <- sapply(
    X = .gaps, FUN = function(.gaps) {
      paste0("\\1", paste0(rep(".", .gaps), collapse = ""), "\\2")
    }
  )

  .cdr3_lengths <- 14:2

  .n1 <- sapply(
    X = .cdr3_lengths,
    FUN = function(.cdr3_lengths) {
      if (.cdr3_lengths %% 2 != 0) (.cdr3_lengths + 1) / 2 - 1 else .cdr3_lengths / 2 - 1
    }
  )

  .n2 <- sapply(
    X = .cdr3_lengths,
    FUN = function(.cdr3_lengths) {
      if (.cdr3_lengths %% 2 != 0) (.cdr3_lengths + 1) / 2 - 2 else .cdr3_lengths / 2 - 1
    }
  )

  .pattern <- mapply(
    FUN = function(.n1, .n2) {
      paste0("(^C.{", .n1, "})(.{", .n2, "}F$|.{", .n2, "}W$)")
    },
    .n1 = .n1,
    .n2 = .n2
  )

  .replace_expressions <- stats::setNames(.replacement, .pattern)

  .seqs <- stringr::str_replace_all(
    .seqs,
    stringr::regex(.replace_expressions, ignore_case = TRUE)
  )

  return(.seqs)
}
