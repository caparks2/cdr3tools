#' Fetch IMGT reference sequences
#'
#' Fetch IMGT reference sequences in FASTA format and convert to a named
#'   character vector.
#'
#' Separate sequence lines under each FASTA entry are concatenated together into
#'   single complete sequences, one for each FASTA entry, with a name taken from
#'   the entry's header.
#'
#' @param .file A character vector. File path or URL to a FASTA formatted
#'   sequence file containing IMGT reference sequences.
#' @param .simplify_headers Logical. `TRUE` (the default) will reduce sequence
#'   names from the full length IMGT FASTA header to the gene name. `FALSE`
#'   leaves the sequence name as the full length IMGT FASTA header.
#' @param .include_species Logical. `TRUE` (the default) will include the
#'   species name in the sequence name. `FALSE` leaves the species name out.
#'
#' @returns A named character vector with length equal to the number of entries
#'   in the input FASTA file.
#'
#' @examplesIf url_exists()
#' # Fetch from IMGT and simplify headers
#' file_path <- paste0(
#'   "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST",
#'   "_reference_directory/Homo_sapiens/TR/TRBV.fasta"
#' )
#'
#' imgt_get_ref_seqs(file_path)
#'
#' @author Christopher Parks
#'
#' @references
#' https://www.imgt.org/vquest/refseqh.html
#'
#' @family IMGT
#' @seealso {imgt_add_junction_gaps()} {imgt_simple_headers()}
#'
#' @export
imgt_get_ref_seqs <- function(.file, .simplify_headers = TRUE,
                              .include_species = .simplify_headers) {

  .fasta <- readLines(.file, warn = FALSE)
  .header_lines <- grep(">", .fasta)

  if (!.simplify_headers) {
    .headers <- .fasta[.header_lines]
  } else {
    .headers <- cdr3tools::imgt_simple_headers(.fasta[.header_lines],
                                               .include_species = .include_species)
  }

  .start <- .header_lines + 1
  .end <- c(.header_lines[2:length(.header_lines)] - 1, length(.fasta))

  .seqs <- mapply(
    function(.start, .end) {
      paste0(.fasta[.start:.end], collapse = "")
    },
    .start = .start,
    .end = .end
  )

  stats::setNames(.seqs, .headers)
}
