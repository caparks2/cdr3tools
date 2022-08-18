#' Format Imunoseq Data for use with Immunarch
#'
#' Prepare Immunoseq data according to Immunarch data format
#'
#' Works on Immunoseq data as returned by the [cdr3tools::read_immunoseq]
#'   function. If unique sequences were collapsed using
#'   [cdr3tools::collapse_sequences] the resulting column `sequence` will be a
#'   list column containing the rearrangment sequences from each of the
#'   duplicate unique sequences.
#' @param .data A list or a single data frame as returned by
#'   [cdr3tools::read_immunoseq].
#' @returns A list or a single data frame of Immunarch formatted data.
#' @author Christopher Parks
#' @references \url{https://immunarch.com/reference/immunr_data_format.html}
#' @examples
#' format_immunarch(cfselo_seqs)
#' @export
format_immunarch <- function(.data) {
  if (!inherits(.data, "list") && inherits(.data, "data.frame")) {
    .data <- list(.data)
  }
  # Format immunoseq files to immunarch format
  immdata <- purrr::map(.data, ~ .x %>%
    # filter for in-frame, productive sequences
    # filter(frame_type == 'In') %>%
    # define CDR3 positions the immunarch way
    dplyr::mutate(
      V.end = n1_index - v_index,
      D.start = d_index - v_index,
      D.end = n2_index - v_index,
      J.start = j_index - v_index
    ) %>%
    # where a region isn't defined, replace the value with NA (validated this with IMGT V-Quest)
    dplyr::mutate(
      dplyr::across(
        .cols = c(V.end, D.start, D.end, J.start),
        .fns = ~ dplyr::case_when(
          .x < 0 ~ NA_real_,
          TRUE ~ .x
        )
      )
    ) %>%
    # select the columns needed for immunarch format. Rename then as we go.
    dplyr::select(
      Clones = templates,
      CDR3.nt = cdr3_rearrangement,
      CDR3.aa = cdr3_amino_acid,
      V.name = v_resolved,
      D.name = d_resolved,
      J.name = j_resolved,
      V.end,
      D.start,
      D.end,
      J.start,
      VD.ins = n1_insertions,
      DJ.ins = n2_insertions,
      Sequence = rearrangement
    ) %>%
    # calculate proportion
    dplyr::mutate(Proportion = Clones / sum(Clones), .after = Clones) %>%
    # set V-J insertions to NA by definition since these are TRB sequences
    dplyr::mutate(VJ.ins = NA, .after = J.start) %>%
    # convert immunoseq VDJ gene names to IMGT gene names
    dplyr::mutate(V.name = cdr3tools::imgt_format_gene_names(V.name)) %>%
    dplyr::mutate(D.name = cdr3tools::imgt_format_gene_names(D.name)) %>%
    dplyr::mutate(J.name = cdr3tools::imgt_format_gene_names(J.name)) %>%
    # sort the data according to abundance
    dplyr::arrange(desc(Proportion)))

  if (length(immdata) == 1) {
    return(immdata[[1]])
  } else {
    return(immdata)
  }
}
