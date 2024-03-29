#' Format TCR VDJ Gene Names with IMGT Rules
#'
#' Reformats character vectors of TCR VDJ gene names to conform to IMGT rules.
#'
#' Vectorised over the input gene names character vector using
#'   stringr::str_replace_all
#'
#' @param .x A character vector of TCR V, D, or J gene names.
#' @param .alleles A logical. `FALSE` (the default) indicates that allele
#'   identifiers (trailing "*01") should be removed. `TRUE` keeps them.
#' @returns A character vector of same length as `.x`.
#' @author Christopher Parks
#' @references
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
#' @examples
#' imgt_format_gene_names(c("TCRBV02-01*01", "TCRBV04-03*01", "TCRBV30-01*01"))
#' @family IMGT
#' @export
imgt_format_gene_names <- function(.x, .alleles = FALSE) {
  res <- stringr::str_replace_all({{ .x }}, ",", ", ")
  res <- stringr::str_replace_all(res, "-([0])([0-9])", "-\\2")
  res <- stringr::str_replace_all(res, "([VDJ])([0])([0-9])", "\\1\\3")
  res <- stringr::str_replace_all(res, "TCR", "TR")
  if (!.alleles) {
    res <- stringr::str_replace_all(res, "\\*[0-9][0-9]", "")
    return(res)
  }
  return(res)
}
