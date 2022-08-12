#' Jensen Shannon Divergence (or Distance)
#'
#' Calculate Jensen Shannon Divergence or Distance between two
#'   antigen receptor sequence repertoires.
#'
#' The Jensen Shannon Divergence (or Distance), JSD, is calculated for normalized
#'   repertoire vectors of counts (or frequencies) using relative entropy
#'   (the Kullback Leibler Divergence). Returned values are between 0 and 1,
#'   with 0 meaning identical repertoires and 1 meaning complete divergence.
#'   JSD is an information theory-based measure of the divergence of TCR
#'   repertoires. This is a symmetric value defined for any two repertoires p
#'   and q as JSD(p, q) = (D(p || m) + D(q || m)) / 2. The code was adapted
#'   as closely as possible from the Python Scipy versions given in the
#'   references.
#'
#' @param p Numeric vector of repertoire counts or frequencies (normalized
#'   counts). Corresponds to the "templates", "Counts", "frequency", or
#'   "Proportion" columns in data returned by [cdr3tools::read_immunoseq] or
#'   [cdr3tools::format_immunarch].
#' @param q Numeric vector of repertoire counts or frequencies (normalized
#'   counts). Corresponds to the "templates", "Counts", "frequency", or
#'   "Proportion" columns in data returned by [cdr3tools::read_immunoseq] or
#'   [cdr3tools::format_immunarch].
#' @param base Numeric. The logarithm base to use in the calculation. Default
#'   is the natural logarithm base e. Note that log base 2 is commonly used.
#' @param distance Logical. `FALSE` (the default) calculates Jensen Shannon
#'   Divergence. `TRUE` calculates Jensen Shannon Distance, the square root of
#'   Jensen Shannon Divergence.
#' @returns Numeric. The value of Jensen Shannon Divergence, or at user's
#'   option, the value of Jensen Shannon Distance, both using the default
#'   natural logarithm or user supplied logarithm base value.
#' @examples
#' p = c(0.10, 0.40, 0.50)
#' q = c(0.80, 0.15, 0.05)
#' # Jensen Shannon Divergence using ln
#' jensen_shannon(p, q)
#' # Jensen Shannon Divergence using log2
#' jensen_shannon(p, q, base = 2)
#' # Jensen Shannon Distance using log2
#' jensen_shannon(p, q, base = 2, distance = TRUE)
#' @author Christopher Parks
#' @references
#' DeWolf S, Grinshpun B, Savage T, Lau SP, Obradovic A, Shonts B, Yang S,
#'   Morris H, Zuber J, Winchester R, Sykes M, Shen Y. Quantifying size and
#'   diversity of the human T cell alloresponse. JCI Insight.
#'   2018 Aug 9;3(15):e121256. doi: 10.1172/jci.insight.121256.
#'   PMID: 30089728; PMCID: PMC6129121.
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6129121/}
#' \url{https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.jensenshannon.html}
#' \url{https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html}
#' \url{https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.rel_entr.html}
#' @export
jensen_shannon <- function(p, q, base = NULL, distance = NULL) {
  if (!any(inherits(p, "numeric"), inherits(q, "numeric")))
    rlang::abort("p and q must be numeric vectors of counts or frequencies.")

  # process inputs
  if (is.null(distance))
    distance <- FALSE
  if (is.null(base))
    base <- exp(1)
  p <- p[p > 0]
  q <- q[q > 0]

  # normalize p and q
  if (!min(p) >= 0 || max(p) <= 1)
    p <- p / sum(p)
  if (!min(q) >= 0 || max(q) <= 1)
    q <- q / sum(q)

  # calculate point-wise mean
  m <- (p + q) / 2

  # calculate relative entropy (Kullback Leibler Divergence)
  Dp <- sum(kullback_leibler(p, m, base))
  Dq <- sum(kullback_leibler(q, m, base))

  # calculate Jensen Shannon Divergence
  jsd <- 0.5 * Dp + 0.5 * Dq

  # calculate Jensen Shannon Distance
  if (distance)
    jsd <- sqrt(jsd)

  return(jsd)
}

kullback_leibler <- function(p, q, base = NULL) {
  # process inputs
  if (is.null(base))
    base <- exp(1)

  # normalize p and q
  p <- p[p > 0]
  q <- q[q > 0]
  if (!min(p) >= 0 || max(p) <= 1)
    p <- p / sum(p)
  if (!min(q) >= 0 || max(q) <= 1)
    q <- q / sum(q)

  # calculate relative entropy (Kullback Leibler Divergence)
  p * log(p / q, base)
}
