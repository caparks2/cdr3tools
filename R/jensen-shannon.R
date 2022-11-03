#' Jensen Shannon Divergence and Distance
#'
#' Calculate Jensen Shannon Divergence or Jensen Shannon Distance between two
#'   antigen receptor sequence repertoires.
#'
#' The Jensen Shannon Divergence (JSD) is calculated for normalized
#'   vectors of unique sequence counts, or frequencies. Returned values are
#'   between 0 and 1, with 0 meaning identical repertoires and 1 meaning
#'   complete divergence. JSD is an information theory measure of the divergence
#'   of TCR repertoires. This is a symmetric value defined for any two
#'   vectors of normalised unique sequence counts *P* and *Q* in:
#'
#'   \deqn{\mathrm{JSD}(P \parallel Q) = \frac{1}{2}D(P \parallel M) + \frac{1}{2}D(Q \parallel M)}
#'
#'   where D is the Kullback Leibler Divergence (also called the relative
#'   entropy) and M is a vector of point-wise means between *P* and *Q*. The
#'   code was adapted as closely as possible from the scipy versions given in
#'   the references.
#'
#'   To reproduce the JSD method used by Boris Grinshpun (Yufeng Shen and Megan
#'   Sykes Labs) in the DeWolf 2018 JCI Insight paper (see references),
#'   set `distance` = `TRUE`.
#'
#'   *A special note:* R recycles atomic vectors during arithmetic operations if
#'   the vectors are of unequal lengths. This means JSD will be calculated
#'   incorrectly in R if the two repertoires being compared have differing
#'   numbers of unique sequences. JSD for antigen receptor repertoires should
#'   not be calculated this way, so the function will exit with an error
#'   message. To prevent this, make sure input vectors p and q are of equal
#'   lengths, and crucially, that each pairwise element of p and q correspond
#'   to strictly the same unique sequence.
#'
#'   `p` and `q` input vectors can be integer counts or normalized count doubles
#'   (frequencies).
#'
#' @param p Numeric. Vector of unique sequence counts or frequencies (normalized
#'   counts). Corresponds to the "templates", "Counts", "frequency", or
#'   "Proportion" columns in data returned by [cdr3tools::read_immunoseq] or
#'   [cdr3tools::format_immunarch].
#' @param q Numeric. Vector of unique sequence counts or frequencies (normalized
#'   counts). Corresponds to the "templates", "Counts", "frequency", or
#'   "Proportion" columns in data returned by [cdr3tools::read_immunoseq] or
#'   [cdr3tools::format_immunarch].
#' @param base Numeric. The logarithm base to use in the calculation. If `NULL`,
#'   the default is base 2. Note that base e and base 10 are also commonly
#'   used.
#' @param distance Logical. If `NULL` or `FALSE` (the default), calculates
#'   Jensen Shannon Divergence. `TRUE` calculates Jensen Shannon Distance, the
#'   square root of Jensen Shannon Divergence.
#' @returns Numeric. The value of Jensen Shannon Divergence, or at user's
#'   option, the value of Jensen Shannon Distance, both using the default
#'   base 2 logarithm or user supplied logarithm base value.
#' @examples
#' p = c(0.10, 0.40, 0.50)
#' q = c(0.80, 0.15, 0.05)
#'
#' # Jensen Shannon Divergence using ln
#' jensen_shannon(p, q, base = exp(1))
#'
#' # Jensen Shannon Divergence using log2
#' jensen_shannon(p, q)
#'
#' # Jensen Shannon Distance using log2
#' jensen_shannon(p, q, base = 2, distance = TRUE)
#' @author Christopher Parks
#' @references
#' DeWolf S, Grinshpun B, Savage T, Lau SP, Obradovic A, Shonts B, Yang S,
#'   Morris H, Zuber J, Winchester R, Sykes M, Shen Y. Quantifying size and
#'   diversity of the human T cell alloresponse. JCI Insight.
#'   2018 Aug 9;3(15):e121256. doi: 10.1172/jci.insight.121256.
#'   PMID: 30089728; PMCID: PMC6129121.
#'
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6129121/}
#'
#' \url{https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.jensenshannon.html}
#'
#' \url{https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.entropy.html}
#'
#' \url{https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.rel_entr.html}
#'
#' @export
jensen_shannon <- function(p, q, base = NULL, distance = NULL) {
  # check for correct input format
  if (!any(inherits(p, "numeric"), inherits(q, "numeric")))
    rlang::abort("p and q must be numeric vectors of counts or frequencies.")
  # process inputs
  if (is.null(distance))  distance  <- FALSE
  if (is.null(base))      base      <- 2
  p                                 <- p[!is.na(p)]
  q                                 <- q[!is.na(q)]
  # normalize p and q if p and q are integer vectors
  if (any(p < 0, q < 0)) rlang::abort("p and q must be positive integers or doubles")
  if (all(p %% 1 == 0)) p <- p / sum(p)
  if (all(q %% 1 == 0)) q <- q / sum(q)
  if (!(all.equal(sum(q), 1) && all.equal(sum(p), 1)))
    rlang::abort(
      paste(
        "p and q must each sum to 1 for i in `p[i] / sum(p)` and",
        "`q[i] / sum(q)`.\n  Please check that `sum(p) == 1` and `sum(q) == 1`",
        "both result in `TRUE`\n  if p and q are vectors of frequencies.",
        "If p and q are vectors of counts\n  check that `sum(p / sum(p)) == 1`",
        "and `sum(q / sum(q)) == 1` both result\n  in `TRUE`."
      )
    )
  # calculate vector of point-wise means (and avoid R's recycling behavior if
  #   input vector lengths are unequal)
  if (length(p) != length(q))
    rlang::abort(
      paste(
        "p and q must have equal lengths, and each pair of p and q must",
        "correspond to the same unique TCR sequence."
      )
    )
  m <- (p + q) / 2
  # calculate JS Divergence
  Dp  <- kullback_leibler(p, m, base)
  Dq  <- kullback_leibler(q, m, base)
  JSD <- 0.5 * (Dp + Dq)
  ## calculate JS Divergence (alternative steps but same answer)
  ## Hj  <- shannons_entropy(m, base)
  ## Hp  <- shannons_entropy(p, base)
  ## Hq  <- shannons_entropy(q, base)
  ## JSD <- Hj - 0.5 * (Hp + Hq)
  # calculate JS Distance or return JS Divergence
  if (distance) sqrt(JSD) else JSD
}

kullback_leibler <- function(p, q, base) {
  x <- p * log(p / q, base)
  x[is.nan(x)] <- 0
  sum(x)
}
