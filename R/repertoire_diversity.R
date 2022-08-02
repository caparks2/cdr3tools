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
#'   If a `.x` is a data frame of list of data frames the column containing
#'   counts must be named one of the following:
#'   \describe{
#'   \item{Column containing counts}{named "reads", "copies", "templates",
#'     "clones", or "counts".}
#'   }
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
#'   See ?diversity_methods
#'
#'   Default is to return all diversity measures.
#'
#' @param .r A numeric vector of length 1. The fraction of top unique sequences
#'   that, in their sum, account for the given `.r` proportion of total copies
#'   (or reads or templates). Only used when `.method = "r20"`. Default value is
#'   0.2.
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
#' @seealso [cdr3tools::diversity_methods()]
#' @family Repertoire Diversity
#' @export
repertoire_diversity <- function(.x, .method = diversity_methods(), .r = 0.20) {

  x <- .x
  method <- .method
  r <- .r

  if (!(inherits(x, "list") || inherits(x, "data.frame") || inherits(x, "numeric"))) {
    rlang::abort(
      paste0(
        "x must either be a single data frame containing a column of sequence ",
        "counts (named `templates`, `reads`, `copies`, `clones`, or `counts`)",
        ", a list containing multiple of such data frames, a single numeric ",
        "vector of sequence counts, or a list containing multiple of such numeric ",
        "vectors."
      )
    )
  }

  if (inherits(x, "list")) {
    result <- lapply(x, function(x) repertoire_diversity_internal(x, method, r))
    result <- do.call(rbind, result)
    if (!is.null(names(x))) {
      result <- cbind(ID = names(x), result)
    }
    result <- tibble::as_tibble(result)
    # names(result) <- method
    return(result)
  }

  if (inherits(x, "data.frame")) {
    result <- repertoire_diversity_internal(x, method, r)
    result <- tibble::as_tibble(result)
    return(result)
  }

  if (inherits(x, "numeric")) {
    result <- repertoire_diversity_internal(x, method, r)
    if (length(names(result)) > 1) {
      result <- tibble::as_tibble(result)
      return(result)
    }
    result <- as.numeric(result)
    return(result)
  }

}

repertoire_diversity_internal <- function(x, method, r) {

  if (inherits(x, "data.frame")) {
    x <- as.data.frame(x)

    if (!any(grepl("templates|reads|copies|clones|counts", names(x), ignore.case = TRUE), na.rm = TRUE)) {
      rlang::abort(
        paste0(
          "x must either be a single data frame containing a column of sequence ",
          "counts (named `templates`, `reads`, `copies`, `clones`, or `counts`)",
          ", a list containing multiple of such data frames, a single numeric ",
          "vector of sequence counts, or a list containing multiple of such numeric ",
          "vectors."
        )
      )
    }

    counts_col <- grep("templates|reads|copies|clones|counts", names(x), ignore.case = TRUE)

    if (length(counts_col) > 1) {
      rlang::abort("There must be only one column containing counts.")
    }

    x <- x[[counts_col]]
  }

  method <- rlang::arg_match(method, diversity_methods(), multiple = TRUE)

  result <- lapply(
    method,
    # FUN.VALUE = numeric(1),
    # USE.NAMES = TRUE,
    FUN = function(method) {
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
        rlang::abort("No valid method detected. see ?diversity_methods")
      )
    }
  )

  names(result) <- method
  result <- as.data.frame(result, row.names = NULL)
  return(result)
}

shannons_entropy <- function(x) {
  x <- x[x > 0]
  S <- length(x)
  P <- x / sum(x)
  Hmax <- log2(S)
  H <- -1 * sum(P * log2(P))
  return(H)
}

shannons_diversity_index <- function(x) {
  x <- x[x > 0]
  S <- length(x)
  P <- x / sum(x)
  Hmax <- log(S)
  H <- -1 * sum(P * log(P))
  return(H)
}

shannons_clonality <- function(x) {
  x <- x[x > 0]
  S <- length(x)
  P <- x / sum(x)
  Hmax <- log2(S)
  H <- -1 * sum(P * log2(P))
  E <- H / Hmax
  C <- 1 - E
  return(C)
}

simpsons_index <- function(x) {
  x <- x[x > 0]
  P <- x / sum(x)
  D1 <- sum(P * P)
  return(D1)
}

gini_simpson_index <- function(x) {
  x <- x[x > 0]
  P <- x / sum(x)
  D1 <- sum(P * P)
  GS <- 1 - D1
  return(GS)
}

simpsons_dominance <- function(x) {
  x <- x[x > 0]
  P <- x / sum(x)
  D1 <- sum(P * P)
  D2 <- 1 / D1
  return(D2)
}

simpsons_equitability <- function(x) {
  x <- x[x > 0]
  S <- length(x)
  P <- x / sum(x)
  D1 <- sum(P * P)
  D2 <- 1 / D1
  E <- D2 / S
  return(E)
}

simpsons_clonality <- function(x) {
  x <- x[x > 0]
  P <- x / sum(x)
  D1 <- sum(P * P)
  return(sqrt(D1))
}

r20 <- function(x, r) {
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
  data <- tibble::tibble(templates = {{ x }}) %>%
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
