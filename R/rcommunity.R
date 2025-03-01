#' Random communities
#'
#' `rcommunity()` draws random communities according to a probability distribution.
#' `rspcommunity()` extends it by spatializing the random communities.
#'
#' Communities of fixed `size` are drawn in a multinomial distribution according
#' to the distribution of probabilities provided by `prob`.
#' An abundance vector `abd` may be used instead of probabilities,
#' then `size` is by default the total number of individuals in the vector.
#' Random communities can be built by drawing in a multinomial law following
#' \insertCite{Marcon2012a;textual}{divent}, or trying to estimate the
#' distribution of the actual community with [probabilities].
#' If `bootstrap` is "Chao2013", the distribution is estimated by a single
#' parameter model and unobserved species are given equal probabilities.
#' If `bootstrap` is "Chao2015", a two-parameter model is used and unobserved
#' species follow a geometric distribution.
#'
#' Alternatively, the probabilities may be drawn following a classical
#' distribution: either lognormal ("lnorm")  \insertCite{Preston1948}{divent}
#' with given standard deviation (`sd_lnorm`; note that the mean is actually
#' a normalizing constant. Its value is set equal to 0 for the simulation of
#' the normal distribution of unnormalized log-abundances), log-series ("lseries")
#' \insertCite{Fisher1943}{divent} with parameter `fisher_alpha`, geometric
#' ("geom") \insertCite{Motomura1932}{divent} with parameter `prob_geom`,
#' or broken stick ("bstick") \insertCite{MacArthur1957}{divent}.
#' The number of simulated species is fixed by `species_number`, except for
#' "lseries" where it is obtained from `fisher_alpha` and `size`:
#' \eqn{S = \alpha \ln(1 + size / \alpha)}.
#' Note that the probabilities are drawn once only.
#' If the number of communities to draw, `n`, is greater than 1, then they are
#' drawn in a multinomial distribution following these probabilities.
#'
#' Log-normal, log-series and broken-stick distributions are stochastic.
#' The geometric distribution is completely determined by its parameters.
#'
#' Spatialized communities include the location of individuals in a window,
#' in a [dbmss::wmppp] object.
#' Several point processes are available, namely binomial (points are uniformly
#' distributed in the window) and \insertCite{Thomas1949;textual}{divent}, which
#' is clustered.
#'
#' Point weights, that may be for instance the size of the trees in a forest
#' community, can be uniform, follow a Weibull or a negative exponential distribution.
#' The latter describe well the diameter distribution of trees in a forest
#' \insertCite{Rennolls1985,Turner2004}{divent}.
#'
#' @inheritParams check_divent_args
#' @param n the number of communities to draw.
#'
#' @name rcommunity
#' @references
#' \insertAllCited{}
NULL


#' @rdname rcommunity
#'
#' @returns `rcommunity()` returns an object of class [abundances].
#' @export
#'
#' @examples
#' # Generate a community made of 100000 individuals among 300 species and fit it
#' abundances <- rcommunity(n = 1, size = 1E5,
#'   species_number = 300, distribution = "lnorm")
#' autoplot(abundances)
rcommunity <- function(
    n,
    size = sum(abd),
    prob = NULL,
    abd = NULL,
    bootstrap = c("Chao2015", "Marcon2012", "Chao2013"),
    species_number = 300,
    distribution = c("lnorm", "lseries", "geom", "bstick"),
    sd_lnorm = 1,
    prob_geom = 0.1,
    fisher_alpha = 40,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    check_arguments = TRUE) {

  # Check arguments
  bootstrap <- match.arg(bootstrap)
  distribution <- match.arg(distribution)
  coverage_estimator <- match.arg(coverage_estimator)
  if (any(check_arguments)) {
    check_divent_args()
    if (!is.null(prob) & !is.null(abd)) {
      cli::cli_abort("{.code prob} and {.code abd} can't be both given.")
    }
  }

  the_prob <- the_abd <- NULL

  # Generate a distribution (prob and abd are null) or use prob or build probabilities from abd
  if (is.null(prob) & is.null(abd)) {
    # Draw in a distribution
    name <- distribution
    # Other distributions: draw probabilities
    the_prob <- switch(distribution,
      geom = prob_geom /
        (1 - (1 - prob_geom) ^ species_number) * (1 - prob_geom) ^ (0:(species_number - 1)),
      lnorm = (the_abd <- stats::rlnorm(species_number, meanlog = 0, sdlog = sd_lnorm)) /
        sum(the_abd),
      lseries = (
        the_abd <-
          rlseries(
            species_number <- fisher_alpha * log(1 + size / fisher_alpha),
            fisher_alpha = fisher_alpha,
            size = size
          )
        ) / sum(the_abd),
      bstick = c(cuts <- sort(stats::runif(species_number - 1)), 1) - c(0, cuts)
    )
  } else {
    # Simulation from a distribution
    name <- ifelse(is.null(prob), "prob", "abd")
  }
  if (!is.null(prob)) {
    # Probabilities are given
    the_prob <- prob[prob != 0]
  }
  if (!is.null(abd)) {
    # Subsample in given abundances.
    # Generate probabilities according to the chosen method.
    if (bootstrap == "Chao2015") {
      the_prob <- probabilities.numeric(
        abd,
        estimator = "Chao2015",
        unveiling = "geometric",
        coverage_estimator = coverage_estimator,
        as_numeric = TRUE,
        check_arguments = FALSE)
    }
    if (bootstrap == "Chao2013") {
      the_prob <- probabilities.numeric(
        abd,
        estimator = "Chao2013",
        unveiling = "uniform",
        coverage_estimator = coverage_estimator,
        as_numeric = TRUE,
        check_arguments = FALSE)
    }
    if (bootstrap == "Marcon2012") {
      the_prob <- abd / sum(abd)
    }
  }

  # Generate communities according to probabilities
  # Draw multinomial samples from probabilities except if abundances have already been obtained (e.g.: lseries)
  the_abd <- t(stats::rmultinom(n, size = size, prob = the_prob))

  return(
    as_abundances.numeric(
      the_abd,
      names = paste(
        name,
        1:n,
        sep = "_"
      ),
      round = FALSE,
      check_arguments = FALSE
    )
  )
}


#' @rdname rcommunity
#' @param spatial the spatial distribution of points.
#' May be "Binomial" (a completely random point pattern except for its fixed number of points) or
#' "Thomas" for a clustered point pattern with parameters `scale` and `mu`.
#' @param species_names a vector of characters or of factors containing the possible species names.
#' @param weight_distribution the distribution of point weights.
#' By default, all weight are 1.
#' May be "uniform" for a uniform distribution between `w_min` and `w_max`,
#' "weibull" with parameters `w_min`, `weibull_scale` and `shape` or
#' "exponential" with parameter `w_mean`.
#'
#' @return `rspcommunity()` returns either a spatialized community,
#' which is a [dbmss::wmppp] object , with `PointType`
#' values as species names if `n`=1 or an object of class ppplist
#' (see [spatstat.geom::solist])  if `n`>1.
#' @export
#'
#' @examples
#' X <- rspcommunity(1, size = 30, species_number = 5)
#' autoplot(X)
#'
rspcommunity <- function(
    n,
    size = sum(abd),
    prob = NULL,
    abd = NULL,
    bootstrap = c("Chao2015", "Marcon2012", "Chao2013"),
    species_number = 300,
    distribution = c("lnorm", "lseries", "geom", "bstick"),
    sd_lnorm = 1,
    prob_geom = 0.1,
    fisher_alpha = 40,
    coverage_estimator = c("ZhangHuang", "Chao", "Turing", "Good"),
    spatial = c("Binomial", "Thomas"),
    thomas_scale = 0.2,
    thomas_mu = 10,
    win = spatstat.geom::owin(),
    species_names = NULL,
    weight_distribution = c("Uniform", "Weibull", "Exponential"),
    w_min = 1,
    w_max = 1,
    w_mean = 15,
    weibull_scale = 20,
    weibull_shape = 2,
    check_arguments = TRUE) {

  # Check arguments
  weight_distribution <- match.arg(weight_distribution)
  bootstrap <- match.arg(bootstrap)
  distribution <- match.arg(distribution)
  coverage_estimator <- match.arg(coverage_estimator)
  spatial <- match.arg(spatial)
  if (any(check_arguments)) {
    check_divent_args()
    if (!is.null(prob) & !is.null(abd)) {
      cli::cli_abort("{.code prob} and {.code abd} can't be both given.")
    }
    if (weight_distribution == "Uniform" & (w_min > w_max)) {
      # Check w_max only in uniform distributions
      cli::cli_abort("{.code w_max} must be greater of equal to {.code w_min}.")
    }
  }

  # Species abundances: call rcommunity
  the_communities <- rcommunity(
    n = n,
    size = size,
    prob = prob,
    abd = abd,
    bootstrap = bootstrap,
    species_number = species_number,
    distribution = distribution,
    sd_lnorm = sd_lnorm,
    prob_geom = prob_geom,
    fisher_alpha = fisher_alpha,
    coverage_estimator = coverage_estimator,
    check_arguments = check_arguments
  )

  # Species names
  is_species_column <- !colnames(the_communities) %in% non_species_columns
  if (is.null(species_names)) {
    # Read the species names returned by rcommunity
    species_names <- names(the_communities)[is_species_column]
  } else {
    # Check the names
    species_names <- species_names[!is.na(species_names)]
    species_names <- unique(species_names)
    species_names <- species_names[nchar(species_names) > 0]
    if (length(species_names) < species_number) {
      cli::cli_abort(
        "The species names vector must contain at least {.code species_number} names."
      )
    } else {
      # Sample species names
      species_names <- sample(species_names, size = species_number)
    }
  }
  # Make them factors
  species_names <- as.factor(species_names)

  # marks_weight() draws random weights for all points of a random point pattern
  marks_weight <- function(size) {
    the_weights <- NULL
    if (weight_distribution == "Uniform") {
      the_weights <- stats::runif(size, min = w_min, max = w_max)
    }
    else if (weight_distribution == "Weibull") {
      the_weights <- w_min +
        stats::rweibull(size, shape = weibull_shape, scale = weibull_scale)
    }
    else if (weight_distribution == "Exponential") {
      the_weights <-  w_min +
        stats::rexp(size, rate = 1 / w_mean)
    }
    return(the_weights)
  }

  # Spatial distribution
  if (spatial == "Binomial") {
    rbinomial <- function(i) {
      # Draw a binomial wmppp
      the_ppp <- spatstat.random::runifpoint(
        the_communities$weight[i],
        win = win
      )
      # Create the marks

      spatstat.geom::marks(the_ppp) <- data.frame(
        PointType = rep(
          species_names,
          times = the_communities[i, is_species_column]
        ),
        PointWeight = marks_weight(the_ppp$n)
      )
      # Set the class
      class(the_ppp) <- c("wmppp", class(the_ppp))
      return(the_ppp)
    }
    # Loop to simulate several point processes
    the_wmppp_list <- lapply(1:n, FUN = rbinomial)
  }
  else if (spatial == "Thomas") {
    rthomas <- function(i) {
      # Prepare an empty point pattern
      the_ppp <- NULL
      # Draw each species
      for (s in 1:species_number) {
        the_ppp_s <- spatstat.random::rThomas(
          kappa = as.numeric(
            the_communities[i, is_species_column][s] /
            spatstat.geom::area.owin(win) / thomas_mu
          ),
          scale = thomas_scale,
          mu = thomas_mu,
          win = win
        )
        # Create the marks
        spatstat.geom::marks(the_ppp_s) <- data.frame(
          PointType = rep(species_names[s], times = the_ppp_s$n),
          PointWeight = marks_weight(the_ppp_s$n)
        )
        # Add the species to the point pattern
        if (is.null(the_ppp)) {
          the_ppp <- the_ppp_s
        } else {
          the_ppp <- spatstat.geom::superimpose(the_ppp, the_ppp_s)
        }
      }
      # Set the class
      class(the_ppp) <- c("wmppp", class(the_ppp))
      return(the_ppp)
    }
    # Loop to simulate several point processes
    the_wmppp_list <- lapply(1:n, FUN = rthomas)
  }

  if (n == 1) {
    # Return a wmppp
    return(the_wmppp_list[[1]])
  } else {
    # Return an object of class ppplist
    return(spatstat.geom::as.solist(the_wmppp_list, promote = TRUE))
  }
}
