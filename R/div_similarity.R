#' Similarity-Based Diversity of a Community
#'
#' Estimate the diversity of species from abundance or probability data and a
#' similarity matrix between species.
#' Several estimators are available to deal with incomplete sampling.
#' Bias correction requires the number of individuals.
#'
#' All species of the `species_distribution` must be found in the matrix of
#' `similarities` if it is named.
#' If it is not, its size must equal the number of species.
#' Then, the order of species is assumed to be the same as that of the
#' `species_distribution`.
#'
#' Similarity-Based diversity can't be interpolated of extrapolated as of the
#' state of the art.
#'
#' @inheritParams check_divent_args
#' @param x An object, that may be a numeric vector containing abundances or probabilities,
#' or an object of class [abundances] or [probabilities].
#' If it is a numeric vector, then its length must equal the dimensions of the
#' `similarities` matrix: species are assumed to be in the same order.
#' @param ... Unused.
#'
#' @returns A tibble with the site names, the estimators used and the estimated diversity.
#'
#' @examples
#' # Similarity matrix
#' Z <- fun_similarity(paracou_6_fundist)
#' # Diversity of each community
#' div_similarity(paracou_6_abd, similarities = Z, q = 2)
#' # gamma diversity
#' div_similarity(paracou_6_abd, similarities = Z, q = 2, gamma = TRUE)
#'
#' @name div_similarity
NULL


#' @rdname div_similarity
#'
#' @export
div_similarity <- function(x, similarities, q = 1, ...) {
  UseMethod("div_similarity")
}


#' @rdname div_similarity
#'
#' @param estimator An estimator of asymptotic diversity.
#'
#' @export
div_similarity.numeric <- function(
    x,
    similarities = diag(length(x)),
    q = 1,
    estimator = c("UnveilJ", "Max", "ChaoShen", "MarconZhang",
                  "UnveilC", "UnveiliC", "naive"),
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    jack_alpha  = 0.05,
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    sample_coverage = NULL,
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  coverage_estimator <- match.arg(coverage_estimator)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
    similarities <- checked_matrix(similarities, x)
  }

  the_entropy <- ent_similarity.numeric(
    x,
    similarities = similarities,
    q = q,
    estimator = estimator,
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    jack_alpha  = jack_alpha,
    jack_max = jack_max,
    coverage_estimator = coverage_estimator,
    sample_coverage = sample_coverage,
    as_numeric = FALSE,
    check_arguments = FALSE
  )
  # Calculate diversity
  the_diversity <- dplyr::mutate(
    the_entropy,
    diversity = exp_q(.data$entropy, q = q),
    .keep = "unused"
  )
  # Return the diversity
  if (as_numeric) {
    return(the_diversity$diversity)
  } else {
    return(the_diversity)
  }
}


#' @rdname div_similarity
#'
#' @export
div_similarity.species_distribution <- function(
    x,
    similarities = diag(sum(!colnames(x) %in% non_species_columns)),
    q = 1,
    estimator = c("UnveilJ", "Max", "ChaoShen", "MarconZhang",
                  "UnveilC", "UnveiliC", "naive"),
    probability_estimator = c("Chao2015", "Chao2013", "ChaoShen", "naive"),
    unveiling = c("geometric", "uniform", "none"),
    jack_alpha  = 0.05,
    jack_max = 10,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    gamma = FALSE,
    as_numeric = FALSE,
    ...,
    check_arguments = TRUE) {

  # Check arguments
  estimator <- match.arg(estimator)
  probability_estimator <- match.arg(probability_estimator)
  unveiling <- match.arg(unveiling)
  coverage_estimator <- match.arg(coverage_estimator)
  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
    similarities <- checked_matrix(similarities, x)
  }

  the_entropy <- ent_similarity.species_distribution(
    x,
    similarities = similarities,
    q = q,
    estimator = estimator,
    probability_estimator = probability_estimator,
    unveiling = unveiling,
    jack_alpha  = jack_alpha,
    jack_max = jack_max,
    coverage_estimator = coverage_estimator,
    as_numeric = FALSE,
    check_arguments = FALSE
  )
  # Calculate diversity
  the_diversity <- dplyr::mutate(
    the_entropy,
    diversity = exp_q(.data$entropy, q = q),
    .keep = "unused"
  )

  # Return the diversity
  if (as_numeric) {
    return(the_diversity$diversity)
  } else {
    return(the_diversity)
  }

}
