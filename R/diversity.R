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
#' @param x Input. A numeric vector of counts. Counts may represent reads,
#'   clone copies, or templates, etc of unique sequences (clonotypes or
#'   rearrangments). Normalization to frequencies is not needed since this is
#'   performed by the functions during calculation of each diversity measure.
#'
#' @param method A character string of one of the following options:
#'   * \code{"shannons_entropy"}
#'   * \code{"shannons_diversity"}
#'   * \code{"shannons_clonality"}
#'   * \code{"simpsons_index"}
#'   * \code{"gini_simpson_index"}
#'   * \code{"simpsons_dominance"}
#'   * \code{"simpsons_equitability"}
#'   * \code{"simpsons_clonality"}
#'   * \code{"r20"}
#'   * \code{"slope"}
#'   The default selection is "simpsons_clonality".
#'
#' @param r A numeric vector of length 1. The fraction of top unique sequences
#'   that, in their sum, account for the given `r` proportion of total copies
#'   (or reads or templates).
#'
#' @returns A numeric vector of length 1.
#'
#' @examples
#' template_counts <- c(100, 8, 3, 2, rep(1, times = 1e4))
#' repertoire_diversity(x = template_counts, method = "r20", r = 0.5)
#'
#' @author Christopher Parks (caparks2@gmail.com)
#'   Boris Grinshpun
#'
#' @references
#' DeWolf S, Grinshpun B, Savage T, Lau SP, Obradovic A, Shonts B, Yang S,
#'   Morris H, Zuber J, Winchester R, Sykes M, Shen Y. Quantifying size and
#'   diversity of the human T cell alloresponse. JCI Insight.
#'   2018 Aug 9;3(15):e121256. doi: 10.1172/jci.insight.121256. PMID: 30089728;
#'   PMCID: PMC6129121
#'
#' @export
repertoire_diversity <- function(
    x,
    method = c(
      "shannons_entropy",
      "shannons_diversity",
      "shannons_clonality",
      "simpsons_index",
      "gini_simpson_index",
      "simpsons_dominance",
      "simpsons_equitability",
      "simpsons_clonality",
      "r20",
      "slope"
    ),
    r = 0.2) {
  method <- method[1]

  result <- switch(method,
    shannons_entropy = shannons_entropy(x),
    shannons_diversity = shannons_diversity_index(x),
    shannons_clonality = shannons_clonality(x),
    simpsons_index = simpsons_index(x),
    gini_simpson_index = gini_simpson_index(x),
    simpsons_dominance = simpsons_dominance(x),
    simpsons_equitability = simpsons_equitability(x),
    simpsons_clonality = simpsons_clonality(x),
    r20 = r20(x, r),
    slope = abundance_slope(x),
    rlang::abort("Please enter an allowed method. see `?repertoire_diversity`.")
  )

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
