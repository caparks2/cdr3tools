#' Measures of Diversity
#'
#' Numeric measures of species diversity derived from information theory
#'   and ecology as applied to quantifying diversity of antigen receptor
#'   repertoires.
#'
#' Shannon's entropy is used in it's original form with log2, as that seems to
#'   be the most commonly found version in the literature. The result is termed
#'   "H". The form used more often in ecology is called Shannon's diversity
#'   index and uses the natural log. For natural log the result is termed "H'",
#'   otherwise it it "H" for log2 and log10. Shannon's Clonality is given as
#'   1 - Pielou's index, which is also 1 - Shannon's equitability. Values range
#'   from 0-1 with 1 being completely clonal while 0 is completely even.
#'   Simpson's diversity index D1 is given as the probability that two clones
#'   (template sequences) taken at random represent the same clonotype (unique
#'   sequence or rearrangement). It is not often used because the values are not
#'   intuitive. The Gini-Simpson index, often referred to as 'Simpson's Index',
#'   is more commonly used than Simpson's diversity index D1 because the values
#'   are more intuitive. Simpson's dominance D2, which is the inverse of
#'   Simpson's index, is also often referred to as 'Simpson's Index'. It's use
#'   is also more common than Simpson's diversity index D1. Simpson's
#'   equitability is also given. Simpson's clonality is included as a normalized
#'   measure and is potentially more useful than the others as it is more
#'   robust to differences in sample sizes. It's discovery was made before
#'   Simpson's diversity index D1 in 1945 by the economist Albert O. Hirschman
#'   but is attributed to Simpson.
#'
#' R20 is a semi-parametric measure of diversity that was
#'   invented by Boris Grinshpun while working in Yufeng Shen's lab, in
#'   collaboration with Megan Sykes, at Columbia University. The measure
#'   describes the cumulative frequency of the most dominant unique sequences (
#'   clonotypes or rearrangements), sorted from most to least abundant, that
#'   account for 20% (or any other selected value) of the total clones (reads,
#'   templates, or copies) in the repertoire.
#'
#' The abundance plot slope estimates the absolute value of the exponent in the
#'   equation describing the power-law (Pareto) distribution. T cell receptor
#'   repertoires have been shown to follow this distribution, with the vast
#'   majority of unique sequences being found with only one copy, and extremely
#'   few unique sequences being found with many copies. The slope of the line
#'   resulting from linear regression of the log transformed abundance and
#'   frequency variables quantifies the diversity of the repertoire. A steep
#'   slope (higher absolute value) indicates a more diverse repertoire, while a
#'   shallow slope (lower absolute value) indicates a more clonal repertoire.
#'
#' Normalization to frequencies is not needed
#'   since this is performed by the functions during calculation of each
#'   diversity measure.
#'
#' @param .x A data frame of antigen receptor sequencing data, minimally
#'   containing a column with counts. Counts may represent reads, copies, or
#'   templates etc of unique sequences (also called clonotypes or
#'   rearrangments). Alternatively, `.x` may be a list of such data frames, a
#'   single numeric vector of counts, or a list of numeric vectors of counts.
#'
#' @param .col A character string (not case sensitive) - the name of the column
#'   containing counts if `.x` is a data frame (or a list of data frames).
#'   `.col` accepts one of the following:
#'     * "reads"
#'     * "seq_reads"
#'     * "copies"
#'     * "templates"
#'     * "clones"
#'     * "counts"
#'     * Defaults to "templates" if no user supplied value and `.x` is a data
#'       frame or list of data frames.
#'
#' @param .method A character string of one of the following options:
#'   * \code{"shannons.clonality"}
#'   * \code{"shannons.entropy"}
#'   * \code{"shannons.diversity"}
#'   * \code{"gini.simpson"}
#'   * \code{"simpsons.clonality"}
#'   * \code{"simpsons.index"}
#'   * \code{"simpsons.dominance"}
#'   * \code{"simpsons.equitability"}
#'   * \code{"r20"}
#'   * \code{"slope"}
#'   See ?cdr3tools::diversity_methods
#'
#'   Default is to return all diversity measures if no user supplied value.
#'
#' @param .r A numeric vector of length 1. The fraction of top unique sequences
#'   that, in their sum, account for the given `.r` proportion of total copies
#'   (or reads or templates). Only used when `.method = "r20"`. Default value is
#'   0.2 if no user supplied value.
#'
#' @returns A numeric vector of length 1, or a data frame of results.
#'
#' @examples
#' template_counts <- c(100, 8, 3, 2, rep(1, times = 1e4))
#' repertoire_diversity(.x = template_counts, .method = "r20", .r = 0.5)
#'
#' @author Boris Grinshpun
#' @author Aleksandar Obradovic
#' @author Christopher Parks
#'
#'
#' @references
#' DeWolf S, Grinshpun B, Savage T, Lau SP, Obradovic A, Shonts B, Yang S,
#'   Morris H, Zuber J, Winchester R, Sykes M, Shen Y. Quantifying size and
#'   diversity of the human T cell alloresponse. JCI Insight.
#'   2018 Aug 9;3(15):e121256. doi: 10.1172/jci.insight.121256. PMID: 30089728;
#'   PMCID: PMC6129121
#'
#' Obradovic A, Shen Y, Sykes M, Fu J. Integrated Analysis Toolset for Defining
#'   and Tracking Alloreactive T-cell Clones After Human Solid Organ and
#'   Hematopoietic Stem Cell Transplantation. Softw Impacts. 2021 Nov;10:100142.
#'   doi: 10.1016/j.simpa.2021.100142. Epub 2021 Sep 23. PMID: 35291378; PMCID:
#'   PMC8920412.
#' @family Immunoseq
#' @export
repertoire_diversity <- function(.x, .col = NULL,
                                 .method = cdr3tools::diversity_methods(),
                                 .r = NULL) {

  x <- .x
  col <- .col
  method <- .method
  if (is.null(.r)) {
    .r <- 0.2
  }
  if (.r > 1 || .r < 0) {
    rlang::abort(".r must be a numeric in the set [0, 1]")
  }
  r <- .r

  if (!any(inherits(x, "list"), inherits(x, "data.frame"), inherits(x, "numeric"), inherits(x, "integer"))) {
    rlang::abort(
      paste0(
        ".x must either be a single data frame containing a column of sequence ",
        "counts (named `templates`, `reads`, `seq_reads`, `copies`, `clones`, or `counts`)",
        ", a list containing multiple of such data frames, a single numeric ",
        "vector of sequence counts, or a list containing multiple of such numeric ",
        "vectors."
      )
    )
  }

  if (inherits(x, "list")) {
    result <- lapply(x, function(x) repertoire_diversity_internal(x, col, method, r))
    result <- do.call(rbind, result)
    if (!is.null(names(x))) {
      result <- cbind(ID = names(x), result)
    }
    result <- tibble::as_tibble(result)
    return(result)
  }

  if (inherits(x, "data.frame")) {
    result <- repertoire_diversity_internal(x, col, method, r)
    result <- tibble::as_tibble(result)
    return(result)
  }

  if (inherits(x, "numeric")) {
    result <- repertoire_diversity_internal(x, col, method, r)
    if (length(names(result)) > 1) {
      result <- tibble::as_tibble(result)
      return(result)
    }
    result <- as.numeric(result)
    return(result)
  }

}

repertoire_diversity_internal <- function(x, col, method, r) {

  if (inherits(x, "data.frame")) {
    x <- as.data.frame(x)

    if (nrow(x) == 0) {
      x <- x[NA, ]
    }

    if (!any(grepl("templates|reads|seq_reads|copies|clones|counts", names(x),
                   ignore.case = TRUE), na.rm = TRUE)) {
      rlang::abort(
        paste0(
          ".x must either be a single data frame containing a column of sequence ",
          "counts (named `templates`, `reads`, `seq_reads`, `copies`, `clones`, or `counts`)",
          ", not case sensitive",
          ", a list containing multiple of such data frames, a single numeric ",
          "vector of sequence counts, or a list containing multiple of such numeric ",
          "vectors."
        )
      )
    }

    if (is.null(col)) {
      col <- grep("clones|templates", names(x), ignore.case = TRUE, value = TRUE)
    }
    if (!any(grepl("templates|reads|seq_reads|copies|clones|counts", col,
                  ignore.case = TRUE))) {
      rlang::abort('`.col` must be one of "templates", "reads", "seq_reads", "copies", "clones", or "counts". Not case sensitive.')
    }
    counts_col <- col
    if (length(counts_col) > 1) {
      rlang::abort("There must be only one column containing counts.")
    }

    x <- x[[counts_col]]
  }

  if (length(x) == 0) {
    x <- NA_real_
  }

  method <- rlang::arg_match(method, diversity_methods(), multiple = TRUE)

  result <- lapply(method, function(method) {
    switch(method,
      shannons.entropy = shannons_entropy(x),
      shannons.diversity = shannons_diversity_index(x),
      shannons.clonality = shannons_clonality(x),
      simpsons.index = simpsons_index(x),
      gini.simpson = gini_simpson_index(x),
      simpsons.dominance = simpsons_dominance(x),
      simpsons.equitability = simpsons_equitability(x),
      simpsons.clonality = simpsons_clonality(x),
      r20 = r20(x, r),
      slope = abundance_slope(x),
      rlang::abort("No valid method detected. see ?cdr3tools::diversity_methods")
    )
    }
  )

  if (any(grepl("r20", method, ignore.case = TRUE)) && r != 0.2) {
    method[grepl("r20", method, ignore.case = TRUE)] <- paste0("r", r * 100)
  }

  names(result) <- method
  result <- as.data.frame(result, row.names = NULL)
  return(result)
}

shannons_entropy <- function(x, base = NULL) {
  if (length(x) == 1 && is.na(x))
    return(NA_real_)
  if (is.null(base))
    base <- 2
  x <- x[which(!is.na(x))]
  P <- x / sum(x)
  H <- -1 * sum(P * log(P, base))
  return(H)
}

shannons_diversity_index <- function(x) {
  if (length(x) == 1 && is.na(x))
    return(NA_real_)
  x <- x[which(!is.na(x))]
  P <- x / sum(x)
  H <- -1 * sum(P * log(P))
  return(H)
}

shannons_clonality <- function(x, base = NULL) {
  if (length(x) == 1 && is.na(x))
    return(NA_real_)
  if (is.null(base))
    base <- 2
  x <- x[which(!is.na(x))]
  x <- x[x > 0]
  S <- length(x)
  P <- x / sum(x)
  Hmax <- log(S, base)
  H <- -1 * sum(P * log(P, base))
  E <- H / Hmax
  C <- 1 - E
  return(C)
}

simpsons_index <- function(x) {
  if (length(x) == 1 && is.na(x)) {
    return(NA_real_)
  }
  x <- x[which(!is.na(x))]
  P <- x / sum(x)
  D1 <- sum(P * P)
  return(D1)
}

gini_simpson_index <- function(x) {
  if (length(x) == 1 && is.na(x)) {
    return(NA_real_)
  }
  x <- x[which(!is.na(x))]
  P <- x / sum(x)
  D1 <- sum(P * P)
  GS <- 1 - D1
  return(GS)
}

simpsons_dominance <- function(x) {
  if (length(x) == 1 && is.na(x)) {
    return(NA_real_)
  }
  x <- x[which(!is.na(x))]
  P <- x / sum(x)
  D1 <- sum(P * P)
  D2 <- 1 / D1
  return(D2)
}

simpsons_equitability <- function(x) {
  if (length(x) == 1 && is.na(x)) {
    return(NA_real_)
  }
  x <- x[which(!is.na(x))]
  x <- x[x > 0]
  S <- length(x)
  P <- x / sum(x)
  D1 <- sum(P * P)
  D2 <- 1 / D1
  E <- D2 / S
  return(E)
}

simpsons_clonality <- function(x) {
  if (length(x) == 1 && is.na(x)) {
    return(NA_real_)
  }
  x <- x[which(!is.na(x))]
  P <- x / sum(x)
  D1 <- sum(P * P)
  return(sqrt(D1))
}

r20 <- function(x, r) {
  if (length(x) == 1 && is.na(x)) {
    return(NA_real_)
  }
  x <- x[which(!is.na(x))]
  x <- x[x > 0]
  S <- length(x)
  P <- x / sum(x)
  P <- sort(P, decreasing = TRUE)
  EP <- cumsum(P)
  N <- length(which(EP <= r))
  R <- N / S
  return(R)
}

#' @importFrom rlang .data
abundance_slope <- function(x) {
  if (length(x) == 1 && is.na(x)) {
    return(NA_real_)
  }
  x <- x[which(!is.na(x))]
  data <- tibble::tibble(templates = x) %>%
    dplyr::filter(.data$templates > 0) %>%
    dplyr::mutate(template_fraction = .data$templates / sum(.data$templates)) %>%
    dplyr::group_by(.data$template_fraction, .data$templates) %>%
    dplyr::summarise(unique_seqs = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(fraction_unique_seqs = .data$unique_seqs / sum(.data$unique_seqs)) %>%
    dplyr::arrange(.data$templates) %>%
    dplyr::mutate(difference = c(NA, abs(diff(log10(.data$template_fraction)))))

  single_unique_seqs <- which(data$unique_seqs == 1)

  if (length(single_unique_seqs) <= 1) {
    linear_rows <- 1:nrow(data)
  } else if (data$difference[single_unique_seqs][2] < 1.5) {
    cutoff <- which(data$unique_seqs == 1)[2]
    linear_rows <- 1:nrow(data[1:cutoff, ])
  } else {
    cutoff <- which(data$unique_seqs == 1)[1]
    linear_rows <- 1:nrow(data[1:cutoff, ])
  }

  model <- stats::lm(
    formula = log10(unique_seqs) ~ log10(template_fraction),
    data = data,
    subset = linear_rows
  )

  slope <- abs(stats::coefficients(model)[[2]])
  return(slope)
}

#' Diversity Calculation Methods
#'
#' For use with the [cdr3tools::repertoire_diversity()] function.
#'
#' These are the various methods of diversity calculation that can be used to
#'   quantify diversity of T cell receptor sequence repertoires. These methods
#'   must be passed as individual character strings to the `.method` argument in
#'   the [cdr3tools::repertoire_diversity()] function. Available options to
#'   `.method` are:
#'   * \code{"shannons.clonality"}
#'   * \code{"shannons.entropy"}
#'   * \code{"shannons.diversity"}
#'   * \code{"gini.simpson"}
#'   * \code{"simpsons.clonality"}
#'   * \code{"simpsons.index"}
#'   * \code{"simpsons.dominance"}
#'   * \code{"simpsons.equitability"}
#'   * \code{"r20"}
#'   * \code{"slope"}
#' @examples diversity_methods()
#' @author Christopher Parks
#' @family Immunoseq
#' @export
diversity_methods <- function() {
  c("shannons.clonality", "shannons.entropy", "shannons.diversity",
    "gini.simpson", "simpsons.clonality", "simpsons.index",
    "simpsons.dominance", "simpsons.equitability", "r20", "slope")
}
