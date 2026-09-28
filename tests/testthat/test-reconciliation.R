context("test-reconciliation")

test_that("reconciliation", {
  lung_deaths_agg <- lung_deaths_long %>% 
    aggregate_key(key, value = sum(value))
  expect_equal(n_keys(lung_deaths_agg), 3)
  expect_equal(
    lung_deaths_agg$value[1:72], 
    lung_deaths_long$value[1:72] + lung_deaths_long$value[72 + (1:72)]
  )
  expect_output(
    print(lung_deaths_agg$key),
    "<aggregated>"
  )
  expect_output(
    print(lung_deaths_agg),
    "<aggregated>"
  )
  
  skip_if_not_installed("fable")
  
  fit_agg <- lung_deaths_agg %>% 
    model(snaive = fable::SNAIVE(value))
  
  fc_agg <- fit_agg %>% forecast()
  fc_agg_reconciled <- fit_agg %>% reconcile(snaive = min_trace(snaive)) %>% forecast()
  
  expect_equal(
    mean(fc_agg$value),
    mean(fc_agg_reconciled$value)
  )
  expect_failure(
    expect_equal(
      fc_agg$value,
      fc_agg_reconciled$value
    ) 
  )
  
  fit_agg <- lung_deaths_agg %>% 
    model(ses = fable::ETS(value ~ error("A") + trend("A") + season("A")))
  fc_agg <- fit_agg %>% forecast()
  fc_agg_reconciled <- fit_agg %>% reconcile(ses = min_trace(ses)) %>% forecast()
  expect_equal(
    mean(fc_agg_reconciled$value[48 + (1:24)]),
    mean(fc_agg_reconciled$value[(1:24)] + fc_agg_reconciled$value[24 + (1:24)]),
  )
  expect_failure(
    expect_equal(
      fc_agg$value,
      fc_agg_reconciled$value
    )
  )
  
  fc_agg_reconciled <- fit_agg %>% reconcile(ses = min_trace(ses, method = "wls_var")) %>% forecast()
  expect_equal(
    mean(fc_agg_reconciled$value[48 + (1:24)]),
    mean(fc_agg_reconciled$value[(1:24)] + fc_agg_reconciled$value[24 + (1:24)])
  )
  expect_failure(
    expect_equal(
      fc_agg$value,
      fc_agg_reconciled$value
    )
  )
  
  fc_agg_reconciled <- fit_agg %>% reconcile(ses = min_trace(ses, method = "ols")) %>% forecast()
  expect_equal(
    mean(fc_agg_reconciled$value[48 + (1:24)]),
    mean(fc_agg_reconciled$value[(1:24)] + fc_agg_reconciled$value[24 + (1:24)])
  )
  expect_failure(
    expect_equal(
      fc_agg$value,
      fc_agg_reconciled$value
    )
  )
  
  fc_agg_reconciled <- fit_agg %>% reconcile(ses = min_trace(ses, method = "mint_cov")) %>% forecast()
  expect_equal(
    mean(fc_agg_reconciled$value[48 + (1:24)]),
    mean(fc_agg_reconciled$value[(1:24)] + fc_agg_reconciled$value[24 + (1:24)])
  )
  expect_failure(
    expect_equal(
      fc_agg$value,
      fc_agg_reconciled$value
    )
  )
})

test_that("top_down reconciles multi-level hierarchies", {
  skip_if_not_installed("fable")

  # 4-level hierarchy: region (2) x product (3) x channel (2) = 12 bottom series
  set.seed(7213)
  sim_hts <- tidyr::expand_grid(
    index   = tsibble::yearmonth("2020 Jan") + 0:35,
    region  = c("North", "South"),
    product = c("A", "B", "C"),
    channel = c("Online", "Offline")
  ) |>
    dplyr::mutate(sales = rpois(dplyr::n(), lambda = 50)) |>
    tsibble::as_tsibble(key = c(region, product, channel), index = index) |>
    aggregate_key(region / product / channel, sales = sum(sales))

  fc_tbl <- sim_hts |>
    model(snaive = fable::SNAIVE(sales)) |>
    reconcile(td = top_down(snaive, method = "forecast_proportions")) |>
    forecast(h = 3) |>
    dplyr::filter(.model == "td") |>
    as_tibble()

  total <- fc_tbl |>
    dplyr::filter(is_aggregated(region), is_aggregated(product), is_aggregated(channel)) |>
    dplyr::arrange(index)

  # Total must equal sum of first level below
  region_sum <- fc_tbl |>
    dplyr::filter(!is_aggregated(region), is_aggregated(product), is_aggregated(channel)) |>
    dplyr::summarise(.mean = sum(.mean), .by = index) |>
    dplyr::arrange(index)
  expect_equal(total$.mean, region_sum$.mean)

  # Total must equal sum of bottom level
  btm_sum <- fc_tbl |>
    dplyr::filter(!is_aggregated(channel)) |>
    dplyr::summarise(.mean = sum(.mean), .by = index) |>
    dplyr::arrange(index)
  expect_equal(total$.mean, btm_sum$.mean)
})

test_that("middle_out reconciles multi-level hierarchies", {
  skip_if_not_installed("fable")

  # 4-level hierarchy: category (3) x segment (2) x sku (3) = 18 bottom series
  set.seed(4871)
  sim_hts <- tidyr::expand_grid(
    index    = tsibble::yearmonth("2020 Jan") + 0:35,
    category = c("Food", "Tech", "Apparel"),
    segment  = c("Premium", "Budget"),
    sku      = c("S1", "S2", "S3")
  ) |>
    dplyr::mutate(sales = rpois(dplyr::n(), lambda = 50)) |>
    tsibble::as_tsibble(key = c(category, segment, sku), index = index) |>
    aggregate_key(category / segment / sku, sales = sum(sales))

  fc_tbl <- sim_hts |>
    model(snaive = fable::SNAIVE(sales)) |>
    reconcile(mo = middle_out(snaive, split = 1)) |>
    forecast(h = 3) |>
    dplyr::filter(.model == "mo") |>
    as_tibble()

  # Split-level (category) nodes must equal sum of their bottom-level descendants
  cat_agg <- fc_tbl |>
    dplyr::filter(!is_aggregated(category), is_aggregated(segment), is_aggregated(sku)) |>
    dplyr::arrange(index, category)
  cat_from_btm <- fc_tbl |>
    dplyr::filter(!is_aggregated(sku)) |>
    dplyr::summarise(.mean = sum(.mean), .by = c(index, category)) |>
    dplyr::arrange(index, category)
  expect_equal(cat_agg$.mean, cat_from_btm$.mean)

  # Total must equal sum of split level
  total <- fc_tbl |>
    dplyr::filter(is_aggregated(category), is_aggregated(segment), is_aggregated(sku)) |>
    dplyr::arrange(index)
  total_from_cat <- fc_tbl |>
    dplyr::filter(!is_aggregated(category), is_aggregated(segment), is_aggregated(sku)) |>
    dplyr::summarise(.mean = sum(.mean), .by = index) |>
    dplyr::arrange(index)
  expect_equal(total$.mean, total_from_cat$.mean)
})

# --- MinT-Ridge helpers -------------------------------------------------------

# Fixtures for the MinT-Ridge tests. `lung_deaths_agg` is scoped to the first
# test_that() block above, so the MinT-Ridge tests build their own.
mint_ridge_agg <- as_tsibble(cbind(mdeaths, fdeaths)) %>%
  aggregate_key(key, value = sum(value))
mint_ridge_mdl <- list(ses = 1)   # min_trace() only wraps its input

test_that("mint_ridge_grid is data driven and glmnet shaped", {
  # lambda_max is the largest eigenvalue of the scale-free W, and d_max >= 1
  # because mean(diag(W/mean(diag(W)))) == 1
  set.seed(2468)
  for (i in 1:20) {
    A <- crossprod(matrix(rnorm(9 * 40), 9, 40)) / 40
    A <- A / mean(diag(A))
    g <- mint_ridge_grid(A, Tn = 40, n_series = 9)
    d_max <- max(eigen(A, symmetric = TRUE, only.values = TRUE)[["values"]])
    expect_gte(d_max, 1 - 1e-8)
    # 101 points: an exact 0 (plain MinT) plus 100 log-spaced values
    expect_length(g, 101)
    expect_equal(g[1], 0)
    expect_equal(g[101], d_max)
    expect_equal(g, sort(g))
    # T > n_series -> 1e-4 floor; otherwise 1e-2
    expect_equal(g[2], d_max * 1e-4)
    expect_equal(mint_ridge_grid(A, Tn = 9, n_series = 9)[2], d_max * 1e-2)
    # geometric spacing
    expect_equal(diff(log(g[2:101])), rep(diff(log(g[2:101]))[1], 99))
  }
})

test_that("mint_ridge_folds never leak and always respect initial", {
  for (Tn in c(12, 20, 47, 60)) {
    for (k in c(2L, 3L, 5L, 10L)) {
      for (initial in unique(c(1L, 5L, ceiling(Tn / 2), Tn - k))) {
        if (initial < 1 || initial + k > Tn) next
        for (window in c("expanding", "rolling")) {
          folds <- mint_ridge_folds(Tn, k, window, initial)
          expect_length(folds, k)
          for (fold in folds) {
            # in range, non-empty, and disjoint
            expect_gte(min(fold$train), 1L)
            expect_lte(max(fold$train), Tn)
            expect_gte(min(fold$valid), 1L)
            expect_lte(max(fold$valid), Tn)
            expect_length(intersect(fold$train, fold$valid), 0L)
            # `initial` is the *minimum* training window
            expect_gte(length(fold$train), initial)
          }
          # Every observation is either trained on or validated on
          used <- sort(unique(c(
            unlist(lapply(folds, `[[`, "train")),
            unlist(lapply(folds, `[[`, "valid"))
          )))
          expect_equal(used, seq_len(Tn))
          if (window == "expanding") {
            # expanding folds are nested
            for (j in seq_len(k - 1L)) {
              expect_true(all(folds[[j]]$train %in% folds[[j + 1L]]$train))
            }
          } else {
            # rolling folds are exactly `initial` long (bar the first)
            for (j in 2:k) expect_length(folds[[j]]$train, initial)
          }
        }
      }
    }
  }
})

test_that("min_trace validates its cross-validation arguments", {
  expect_error(min_trace(mint_ridge_mdl, k = 1), "at least 2")
  expect_error(min_trace(mint_ridge_mdl, k = 0), "at least 2")
  expect_error(min_trace(mint_ridge_mdl, k = 2.5), "at least 2")
  expect_error(min_trace(mint_ridge_mdl, k = c(3, 4)), "at least 2")
  expect_error(min_trace(mint_ridge_mdl, window = "nope"))
  expect_error(min_trace(mint_ridge_mdl, grid = c(-1, 1)), "non-negative")
  expect_error(min_trace(mint_ridge_mdl, grid = c(0, NA)), "non-negative")
  expect_error(min_trace(mint_ridge_mdl, grid = numeric(0)), "non-negative")
  expect_error(min_trace(mint_ridge_mdl, grid = "a"), "non-negative")

  # the default method is unchanged
  expect_equal(attr(min_trace(mint_ridge_mdl), "method"), "wls_var")
  expect_equal(attr(min_trace(mint_ridge_mdl, "mint_cov"), "method"), "mint_cov")
  expect_equal(attr(min_trace(mint_ridge_mdl, method = "mint_ridge"), "method"), "mint_ridge")

  # CV arguments are stored, sorted and deduplicated
  expect_equal(
    attr(min_trace(mint_ridge_mdl, method = "mint_ridge", grid = c(3, 1, 2, 1)), "grid"),
    c(1, 2, 3)
  )

  # ... and warn when they cannot be used
  expect_warning(min_trace(mint_ridge_mdl, k = 7), "mint_ridge")
  expect_warning(min_trace(mint_ridge_mdl, window = "rolling"), "mint_ridge")
  expect_warning(min_trace(mint_ridge_mdl, grid = 1), "mint_ridge")
  expect_silent(min_trace(mint_ridge_mdl, method = "mint_cov"))
  expect_silent(min_trace(mint_ridge_mdl, method = "mint_ridge"))
})

# --- MinT-Ridge reconciliation ------------------------------------------------

test_that("min_trace(method = 'mint_ridge') reconciles", {
  skip_if_not_installed("fable")

  fit_agg <- mint_ridge_agg %>%
    model(ses = fable::ETS(value ~ error("A") + trend("A") + season("A")))

  for (window in c("expanding", "rolling")) {
    for (sparse in c(TRUE, FALSE)) {
      fc <- fit_agg %>%
        reconcile(ses = min_trace(ses, method = "mint_ridge", window = window, sparse = sparse)) %>%
        forecast()
      # aggregate stays the sum of the two bottom series
      expect_equal(
        mean(fc$value[48 + (1:24)]),
        mean(fc$value[(1:24)] + fc$value[24 + (1:24)])
      )
      # and reconciliation is not a no-op
      expect_failure(expect_equal(
        fit_agg %>% forecast() %>% pull(value),
        fc %>% pull(value)
      ))
    }
  }
})

test_that("mint_ridge spans mint_cov (alpha=0) to OLS (alpha=inf), not bottom-up", {
  skip_if_not_installed("fable")

  fit_agg <- mint_ridge_agg %>%
    model(ses = fable::ETS(value ~ error("A") + trend("A") + season("A")))

  # alpha = 0 is exactly the unpenalised estimator, mint_cov
  expect_equal(
    fit_agg %>% reconcile(ses = min_trace(ses, method = "mint_cov")) %>% forecast(),
    fit_agg %>% reconcile(ses = min_trace(ses, method = "mint_ridge", grid = 0)) %>% forecast()
  )

  # alpha -> inf is exactly OLS, i.e. W = I
  ols_fc <- fit_agg %>% reconcile(ses = min_trace(ses, method = "ols")) %>% forecast()
  expect_equal(
    ols_fc,
    fit_agg %>% reconcile(ses = min_trace(ses, method = "mint_ridge", grid = 1e10)) %>% forecast()
  )

  # ... and OLS is *not* bottom-up. Bottom-up passes the bottom forecasts
  # through untouched, whereas OLS mixes siblings and can use negative weights,
  # so the large-alpha limit of the ridge path must differ from bottom-up.
  bu_fc <- fit_agg %>% reconcile(ses = bottom_up(ses)) %>% forecast()
  expect_failure(expect_equal(ols_fc, bu_fc))

  # Algebraically, for a 2-level hierarchy S = [I; 1 1]:
  #   P_ols = (S'S)^-1 S' = (1/3)[2 -1; -1 2] S'   vs   P_bu = [I 0]
  S <- rbind(c(1, 0), c(0, 1), c(1, 1))
  P_ols <- solve(t(S) %*% S) %*% t(S)
  P_bu <- cbind(diag(2), matrix(0, 2, 1))   # bottom forecasts passed through unchanged
  expect_false(isTRUE(all.equal(P_ols, P_bu)))
  expect_equal(max(abs(P_ols - P_bu)), 1 / 3)
  # both are coherent, they are just different projections
  expect_equal(P_ols %*% S, diag(2), ignore_attr = TRUE)
  expect_equal(P_bu %*% S, diag(2), ignore_attr = TRUE)
})

test_that("mint_ridge is deterministic, agrees sparse/dense, and is customisable", {
  skip_if_not_installed("fable")

  fit_agg <- mint_ridge_agg %>%
    model(ses = fable::ETS(value ~ error("A") + trend("A") + season("A")))

  f_ridge <- function(...) {
    fit_agg %>%
      reconcile(ses = min_trace(ses, method = "mint_ridge", ...)) %>%
      forecast()
  }

  # fully deterministic: the search uses no randomness
  expect_equal(f_ridge(), f_ridge())
  # the sparse and dense reconciliation paths agree
  expect_equal(f_ridge(sparse = TRUE), f_ridge(sparse = FALSE))
  # fold count and grid size change the penalty that gets selected
  expect_equal(mean(f_ridge(k = 2, grid = c(0, 0.1, 5)) %>% pull(value)),
               mean(f_ridge(k = 2, grid = c(0, 0.1, 5)) %>% pull(value)))
  for (grid in list(c(0, 0.5), c(2, 5, 20), 1e-3)) {
    fc <- f_ridge(grid = grid)
    expect_equal(
      mean(fc$value[48 + (1:24)]),
      mean(fc$value[(1:24)] + fc$value[24 + (1:24)])
    )
  }
})

test_that("mint_ridge works where mint_cov is undefined", {
  skip_if_not_installed("fable")
  skip_if_not_installed("tidyr")

  # More bottom series than observations, so the sample covariance is singular
  set.seed(1357)
  wide <- tidyr::expand_grid(
    index = tsibble::yearmonth("2018 Jan") + 0:59,
    btm  = sprintf("s%03d", 1:100)
  ) |>
    dplyr::mutate(value = rnorm(dplyr::n(), 50, 10)) |>
    tsibble::as_tsibble(key = btm, index = index) |>
    aggregate_key(btm, value = sum(value))

  fit_wide <- wide %>% model(snaive = fable::SNAIVE(value))

  # mint_cov aborts on a non positive definite covariance matrix
  expect_error(
    fit_wide %>% reconcile(snaive = min_trace(snaive, method = "mint_cov")) %>% forecast(),
    "positive definite"
  )

  # mint_ridge regularises instead, warns, and returns a coherent forecast
  expect_warning(
    fc <- fit_wide %>%
      reconcile(snaive = min_trace(snaive, method = "mint_ridge")) %>%
      forecast(h = 2),
    "not positive definite"
  )
  expect_s3_class(fc, "fbl_ts")
  expect_equal(nrow(fc), 101L * 2L)   # 100 bottom series + 1 aggregate
  expect_false(anyNA(fc$.mean))
})

test_that("mint_ridge rejects infeasible cross-validation splits", {
  skip_if_not_installed("fable")

  fit_agg <- mint_ridge_agg %>%
    model(snaive = fable::SNAIVE(value))

  # 72 observations, so k = 100 folds cannot be built
  expect_error(
    fit_agg %>%
      reconcile(snaive = min_trace(snaive, method = "mint_ridge", k = 100)) %>%
      forecast(),
    "complete in-sample observations"
  )
  # `initial` must leave room for k validation observations
  expect_error(
    fit_agg %>%
      reconcile(snaive = min_trace(snaive, method = "mint_ridge", k = 5, initial = 70)) %>%
      forecast(),
    "must be between 1 and"
  )
  expect_error(
    fit_agg %>%
      reconcile(snaive = min_trace(snaive, method = "mint_ridge", initial = 0)) %>%
      forecast(),
    "must be between 1 and"
  )
})
