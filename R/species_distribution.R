#' Species Distributions
#'
#' A Species Distribution is a [tibble::tibble] containing species abundances
#' or probabilities.

#' Rows of the tibble are communities and column are species.
#' Values are either abundances or probabilities.
#' Special columns contain the site names, and their weights
#' (e.g. their area or number of individuals):
#' their names must be "site" and "weight".
#' All other column names are considered as species names.
#'
#' `species_distribution` objects include `abundances` and `probabilities`
#' objects.
#'
#' `species_distribution()` creates a `species_distribution` object from a vector
#' or a matrix or a dataframe.
#'
#' `as_species_distribution()`, `as_abundances()` and `as_probabilities` format
#' the numeric, matrix or dataframe `x` so that appropriate
#' versions of community functions (generic methods such as [plot] or
#' [div_richness]) are applied.
#' Abundance values are rounded (by default) to the nearest integer.
#' They also accept a [dbmss::wmppp] objects,
#' i.e. a weighted, marked planar point pattern and count the abundances of
#' point types, character and factor objects.
#'
#' `as_probabilities()` normalizes the vector `x` so that it sums to 1. It gives
#' the same output as `probabilities()` with `estimator = "naive"`.
#'
#' `species_distribution` objects objects can be plotted by [plot] and [autoplot].
#'
#' @inheritParams check_divent_args
#' @param x an object.
#' @param ... Unused.

#' @returns An object of classes `species_distribution` and `abundances`
#' or `probabilities`.
#'
#' `as.double()` and its synonymous `as.numeric()` return a numeric vector
#' that contains species abundances or probabilities of a single-row
#' `species_distribution`.
#' `as.matrix()` returns a numeric matrix if the `species_distribution` contains
#' several rows.
#' These are methods of the generic functions for class `species_distribution`.
#'
#' @examples
#' # Paracou data is a tibble
#' paracou_6_abd
#' # Class
#' class(paracou_6_abd)
#' is_species_distribution(paracou_6_abd)
#' # Whittaker plot fitted by a log-normal distribution
#' autoplot(paracou_6_abd[1,], fit_rac = TRUE, distribution = "lnorm")
#' # Character vectors
#' as_abundances(c("A", "C", "B", "C"))
#'
#' @name species_distribution
NULL


#  Species Distribution ----

#' @rdname species_distribution
#'
#' @param names The names of the species distributions.
#'
#' @export
species_distribution <- function(
    x,
    names = NULL,
    weights = NULL,
    check_arguments = TRUE) {

  # Check the data ----
  if (any(check_arguments)) {
    check_divent_args()
    if (!is.numeric(x)) {
      cli::cli_abort("{.code x} must be numeric")
    }
    if (length(dim(x)) > 2) {
      cli::cli_abort("{.code x} may be a vector or a matrix")
    }
    if (any(x < 0)) {
      cli::cli_abort("Species probabilities or abundances must be positive.")
    }
  }

  # Build a tibble from the data ----
  if (is.vector(x)) {
    ## Single distribution ----
    if (is.null(names(x))) {
      ### Columns: add default species names such as sp_1 ----
      names(x) <- paste(
        "sp",
        formatC(seq_along(x), width = ceiling(log10(length(x))), flag = "0"),
        sep = "_"
      )
    }
    if (length(names) != 1) {
      ### Rows: Add a site name ----
      names <- "site_1"
    }
    # Build a tibble
    the_distribution <- tibble::tibble_row(
      site = names,
      weight = sum(x),
      as.data.frame(t(x))
    )
  } else {
    ## Several distributions ----
    if (is.null(colnames(x))) {
      ### Columns: add default species names such as sp_1 ----
      colnames(x) <- paste(
        "sp",
        formatC(seq_len(ncol(x)), width = ceiling(log10(ncol(x))), flag = "0"),
        sep = "_"
      )
    }
    # Build a tibble
    the_distribution <- tibble::as_tibble(x, rownames = "site")
    ### Rows: site names = names or matrix row names or default ----
    if (!is.null(names)) {
      # site = names if the size matches
      if (length(names) == nrow(x)) {
        the_distribution$site <- names
      } else {
        cli::cli_abort(
          paste(
            "The length of {.code names} must match",
            "the number of lines of the data matrix."
          )
        )
      }
    } else {
      # names is null...
      if (is.null(row.names(x))) {
        # ...and no row names: set default names such as site_1
        the_distribution$site <- paste(
          "site",
          formatC(seq_len(nrow(x)), width = ceiling(log10(nrow(x))), flag = "0"),
          sep = "_"
        )
      }
    }
    ### Rows: site weights ----
    if (!is.null(weights)) {
      # site = weights if the size matches
      if (length(weights) == nrow(x)) {
        the_distribution <- tibble::add_column(
          the_distribution,
          weight = weights,
          .after = "site"
        )
      } else {
        cli::cli_abort(
          paste(
            "The length of {.code weights} must match",
            "the number of lines of the data matrix."
          )
        )
      }
    } else {
      # Weights are the number of individuals
      the_distribution <- tibble::add_column(
        the_distribution,
        weight = rowSums(x),
        .after = "site"
      )
    }
  }

  # Set the class and return ----
  class(the_distribution) <- c("species_distribution", class(the_distribution))
  return(the_distribution)
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution <- function(x, ...) {
  UseMethod("as_species_distribution")
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.numeric <- function(
    x,
    ...,
    check_arguments = TRUE) {

  return(
    species_distribution(
      x,
      check_arguments = check_arguments
    )
  )
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.matrix <- function(
    x,
    names = NULL,
    weights = NULL,
    ...,
    check_arguments = TRUE) {

  return(
    species_distribution(
      x,
      names = names,
      weights = weights,
      check_arguments = check_arguments
    )
  )
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.data.frame <- function(
    x,
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) check_divent_args()
  # Check the data
  if (any(x < 0)) {
    cli::cli_abort("All numeric values of the dataframe must be positive.")
  }

  # Build a tibble
  the_distribution <- tibble::as_tibble(x)

  # The first column should be "site"
  if (!"site" %in% colnames(the_distribution)) {
    the_distribution <- tibble::add_column(
      the_distribution,
      site = paste(
        "site",
        formatC(
          seq_len(nrow(the_distribution)),
          width = ceiling(log10(nrow(the_distribution))),
          flag = "0"
        ),
        sep = "_"
      ),
      .before = 1
    )
  }

  # The second column should be "weight"
  if (!"weight" %in% colnames(the_distribution)) {
    the_distribution <- tibble::add_column(
      the_distribution,
      weight = rowSums(the_distribution[colnames(the_distribution) != "site"]),
      .after = "site"
    )
  }

  # Set the class and return
  class(the_distribution) <- c("species_distribution", class(the_distribution))
  return(the_distribution)
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.wmppp <- function(
    x,
    ...,
    check_arguments = TRUE) {

  return(
    species_distribution(
      as_named_vector.wmppp(x),
      check_arguments = check_arguments
    )
  )
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.character <- function(
    x,
    ...,
    check_arguments = TRUE) {

  # Count the number of items by type
  return(
    species_distribution(
      as_named_vector.character(x),
      check_arguments = check_arguments
    )
  )
}


#' @rdname species_distribution
#'
#' @export
as_species_distribution.factor <- function(
    x,
    ...,
    check_arguments = TRUE) {

  # Count the number of items by type
  return(
    species_distribution(
      # tapply keep factors' names
      tapply(x, x, length) ,
      check_arguments = check_arguments
    )
  )
}


#' @rdname species_distribution
#'
#' @export
is_species_distribution <- function(x) {
  inherits(x, "species_distribution")
}



#  Probabilities ----

#' @rdname species_distribution
#'
#' @export
as_probabilities <- function(x, ...) {
  UseMethod("as_probabilities")
}


#' @rdname species_distribution
#'
#' @export
as_probabilities.numeric <- function(
    x,
    ...,
    check_arguments = TRUE) {

  if (any(x < 0)) {
    cli::cli_abort("Species probabilities must be positive.")
  }

  # Normalize to 1
  prob <- x / sum(x)
  the_probabilities <- as_species_distribution(
    prob,
    check_arguments = check_arguments
  )

  class(the_probabilities) <- c("probabilities", class(the_probabilities))
  return(the_probabilities)
}


#' @rdname species_distribution
#'
#' @export
as_probabilities.matrix <- function(
    x,
    names = NULL,
    weights = NULL,
    ...,
    check_arguments = TRUE) {

  # Calculate probabilities by row
  prob <- x / rowSums(x)

  # Build the species distribution
  the_probabilities <- as_species_distribution.matrix(
    prob,
    names = names,
    weights = weights,
    check_arguments = check_arguments
  )

  # Set the class
  class(the_probabilities) <- c("probabilities", class(the_probabilities))
  return(the_probabilities)
}


#' @rdname species_distribution
#'
#' @export
as_probabilities.data.frame <- function(
    x,
    ...,
    check_arguments = TRUE) {

  # Build the species distribution to add site and weight columns if needed
  abundances <- as_species_distribution.data.frame(
    x,
    check_arguments = check_arguments
  )

  # Select species columns
  species_columns <- !colnames(abundances) %in% non_species_columns
  # Extract abundances
  abd <- as.matrix(abundances[, species_columns])
  # Normalize them
  prob <- abd / rowSums(abd)
  # Build the tibble
  the_probabilities <- cbind(
    data.frame(site = abundances$site, weight = abundances$weight),
    as.data.frame(prob)
  )
  # Restore exact species names (spaces may have been transformed into "_")
  colnames(the_probabilities[, species_columns]) <- colnames(abundances[, species_columns])

  # Build the species distribution again for the classes
  the_probabilities <- as_species_distribution.data.frame(
    the_probabilities,
    check_arguments = FALSE
  )

  # Set the class
  class(the_probabilities) <- c("probabilities", class(the_probabilities))
  return(the_probabilities)
}


#' @rdname species_distribution
#'
#' @export
as_probabilities.wmppp <- function(
    x,
    ...,
    check_arguments = TRUE) {

  the_probabilities <- as_species_distribution(
    as_named_vector.wmppp(x),
    check_arguments = check_arguments
  )

  class(the_probabilities) <- c("probabilities", class(the_probabilities))
  return(the_probabilities)
}


#' @rdname species_distribution
#'
#' @export
as_probabilities.character <- function(
    x,
    ...,
    check_arguments = TRUE) {

  the_probabilities <- as_species_distribution(
    as_named_vector.character(x),
    check_arguments = check_arguments
  )

  class(the_probabilities) <- c("probabilities", class(the_probabilities))
  return(the_probabilities)
}


#' @rdname species_distribution
#'
#' @export
as_probabilities.factor <- function(
    x,
    ...,
    check_arguments = TRUE) {

  the_probabilities <- as_species_distribution(
    # tapply keep factors' names
    tapply(x, x, length) ,
    check_arguments = check_arguments
  )

  class(the_probabilities) <- c("probabilities", class(the_probabilities))
  return(the_probabilities)
}


#' @rdname species_distribution
#'
#' @export
is_probabilities <- function(x) {
  inherits(x, "probabilities")
}


#  Abundances ----

#' @rdname species_distribution
#'
#' @param round If `TRUE`, the values of `x` are rounded to the nearest integer.
#'
#' @export
abundances <- function(
    x,
    round = TRUE,
    names = NULL,
    weights = NULL,
    check_arguments = TRUE) {

  if (any(check_arguments)) {
    check_divent_args()
    if (!is.numeric(x)) {
      cli::cli_abort("{.code x} must be numeric")
    }
    if (any(x < 0)) {
      cli::cli_abort("Species abundances must be positive.")
    }
  }

  if (round) {
    x <- round(x)
  }

  abundances <- species_distribution(
    x,
    names = names,
    weights = weights,
    check_arguments = FALSE
  )

  class(abundances) <- c("abundances", class(abundances))
  return(abundances)
}

#' @rdname species_distribution
#'
#' @export
as_abundances <- function(x, ...) {
  UseMethod("as_abundances")
}


#' @rdname species_distribution
#'
#' @export
as_abundances.numeric <- function(
    x,
    round = TRUE,
    ...,
    check_arguments = TRUE) {

  if (any(check_arguments)) {
    check_divent_args()
    if (any(x < 0)) {
      cli::cli_abort("Species abundances must be positive.")
    }
  }

  if (round) {
    the_names <- names(x)
    x <- round(x)
    names(x) <- the_names
  }

  the_abundances <- as_species_distribution(x, check_arguments = FALSE)

  class(the_abundances) <- c("abundances", class(the_abundances))
  return(the_abundances)
}


#' @rdname species_distribution
#'
#' @export
as_abundances.matrix <- function(
    x,
    round = TRUE,
    names = NULL,
    weights = NULL,
    ...,
    check_arguments = TRUE) {

  if (round) {
    x <- round(x)
  }

  the_abundances <- as_species_distribution.matrix(
    x,
    names = names,
    weights = weights,
    check_arguments = check_arguments
  )

  class(the_abundances) <- c("abundances", class(the_abundances))
  return(the_abundances)
}


#' @rdname species_distribution
#'
#' @export
as_abundances.data.frame <- function(
    x,
    ...,
    check_arguments = TRUE) {

  the_abundances <- as_species_distribution.data.frame(
    x,
    check_arguments = check_arguments
  )

  class(the_abundances) <- c("abundances", class(the_abundances))
  return(the_abundances)
}


#' @rdname species_distribution
#'
#' @export
as_abundances.wmppp <- function(
    x,
    ...,
    check_arguments = TRUE) {

  the_abundances <- as_species_distribution(
    as_named_vector.wmppp(x),
    check_arguments = check_arguments
  )

  class(the_abundances) <- c("abundances", class(the_abundances))
  return(the_abundances)
}


#' @rdname species_distribution
#'
#' @export
as_abundances.character <- function(
    x,
    ...,
    check_arguments = TRUE) {

  the_abundances <- as_species_distribution(
    as_named_vector.character(x),
    check_arguments = check_arguments
  )

  class(the_abundances) <- c("abundances", class(the_abundances))
  return(the_abundances)
}


#' @rdname species_distribution
#'
#' @export
as_abundances.factor <- function(
    x,
    ...,
    check_arguments = TRUE) {

  the_abundances <- as_species_distribution(
    # tapply keep factors' names
    tapply(x, x, length) ,
    check_arguments = check_arguments
  )

  class(the_abundances) <- c("abundances", class(the_abundances))
  return(the_abundances)
}


#' @rdname species_distribution
#'
#' @export
is_abundances <- function(x) {
  inherits(x, "abundances")
}


#  Opposite conversions ----

#' @rdname species_distribution
#'
#' @param use.names If `TRUE`, the names of the `species_distribution` are kept
#' in the matrix or vector they are converted to.
#'
#' @export
as.matrix.species_distribution <- function(x, use.names = TRUE, ...) {
  # Delete co
  the_matrix <- as.matrix.data.frame(x[, !colnames(x) %in% non_species_columns])
  if (!use.names) {
    rownames(the_matrix) <- colnames(the_matrix) <- NULL
  }
  return(the_matrix)
}


#' @rdname species_distribution
#'
#' @export
as.double.species_distribution <- function(x, use.names = TRUE, ...) {
  if (nrow(x) > 1) {
    cli::cli_alert_warning(
      "The species_distribution object contains several rows."
    )
    cli::cli_alert("{.fn as.matrix} is used.")
    return(as.matrix.species_distribution(x, use.names, ...))
  } else {
    is_species_column <- !colnames(x) %in% non_species_columns
    the_vector <- unlist(x[1, is_species_column], use.names = use.names)
    return(the_vector)
  }
}


#' @rdname species_distribution
#'
#' @export
as.numeric.species_distribution <- function(x, use.names = TRUE, ...) {
  return(as.double.species_distribution(x, use.names, ...))
}


#  Utilities ----

#' as_named_vector.character
#'
#' Counts the number of points of a `character` vector and returns a named vector.
#' Names are the items of the character vector.
#' This is equivalent to `as.numeric(table(x))` but `table()`
#' looses the names.
#'
#' @param x a character vector.
#'
#' @returns A named vector with the number of items by name.
#' @noRd
#'
as_named_vector.character <- function(x){
  # Count the number of items. Returns a 1D array, not a vector.
  the_array <- tapply(x, x, length)
  the_vector <- as.vector(the_array)
  # Add the names
  names(the_vector) <- names(the_array)
  return(the_vector)
}


#' as_named_vector.wmppp
#'
#' Counts the number of points of a `wmppp` object and returns a named vector.
#' Names are the point types.
#' This is equivalent to `as.numeric(table(X$marks$PointType))` but `table()`
#' looses the names.
#'
#' @param X a [dbmss::wmppp] object, i.e. a weighted, marked planar point pattern.
#'
#' @returns A named vector with the number of points by type.
#' @noRd
#'
as_named_vector.wmppp <- function(X){
  # Count the number of points by type
  the_vector <- as_named_vector.character(spatstat.geom::marks(X)$PointType)
  # Eliminate NAs due to factor levels with no item
  return(the_vector[!is.na(the_vector)])
}
