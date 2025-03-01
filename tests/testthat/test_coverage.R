# Check bad arguments
testthat::test_that(
  "Warning is returned", {
    testthat::skip_on_cran()
    # Singletons only
    testthat::expect_message(
      coverage(rep(1,5)),
      "Sample coverage is 0, most estimators will return `NaN`."
    )
  }
)

# Check rarely used estimators
testthat::test_that(
  "Coverage is estimated", {
    testthat::skip_on_cran()
    # Chao
    testthat::expect_lt(
      abs(coverage(seq_len(5), estimator="Chao")$coverage - coverage(seq_len(5))$coverage),
      1/1000
    )
    # Turing
    testthat::expect_lt(
      abs(coverage(seq_len(5), Estimator="turing")$coverage - coverage(seq_len(5))$coverage),
      1/100
    )
  }
)

# Combine all parameters
abundances <- paracou_6_abd[1, ]

testthat::test_that(
  "No estimator fails", {
    testthat::skip_on_cran()
    # Estimate coverage systematically
    coverage.list <- lapply(
      # All estimators
      eval(formals(divent:::coverage.numeric)$estimator),
      function(estimator) {
        # print(estimator)
        suppressWarnings(
          coverage(
            abundances,
            estimator = estimator,
            level = NULL,
            as_numeric = FALSE,
            check_arguments = TRUE
          )
        )
      }
    )
    # Coerce to a dataframe
    coverage.dataframe <- do.call(rbind, coverage.list)
    # All values must be < 1
    testthat::expect_gt(
      1,
      max(coverage.dataframe$coverage)
    )
  }
)

# Interpolation and extrapolation
sample_size <- abd_sum(abundances, as_numeric = TRUE)
levels <- c(sample_size / 2, round(sample_size * 1.5))

testthat::test_that(
  "No estimator fails during interpolation and extrapolation", {
    testthat::skip_on_cran()
    # Estimate coverage systematically
    coverage.list <- lapply(
      # All estimators
      eval(formals(divent:::coverage.numeric)$estimator),
      function(estimator) {
        the_list <-lapply(
          # Two levels
          levels,
          function(level) {
            # print(paste(estimator, level))
            suppressWarnings(
              coverage(
                abundances,
                estimator = estimator,
                level = level,
                as_numeric = FALSE,
                check_arguments = TRUE
              )
            )
          }
        )
        # Make a dataframe with the list to avoid nested lists
        the_df <- do.call(rbind, the_list)
      }
    )
    # Coerce to a dataframe
    coverage.dataframe <- do.call(rbind, coverage.list)
    # All values must be < 1
    testthat::expect_gt(
      1,
      max(coverage.dataframe$coverage)
    )
  }
)
