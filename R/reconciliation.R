#' Forecast reconciliation 
#' 
#' This function allows you to specify the method used to reconcile forecasts
#' in accordance with its key structure.
#' 
#' @param .data A mable.
#' @param ... Reconciliation methods applied to model columns within `.data`.
#' 
#' @examplesIf requireNamespace("fable", quietly = TRUE)
#' library(fable)
#' lung_deaths_agg <- as_tsibble(cbind(mdeaths, fdeaths)) %>%
#'   aggregate_key(key, value = sum(value))
#' 
#' lung_deaths_agg %>%
#'   model(lm = TSLM(value ~ trend() + season())) %>%
#'   reconcile(lm = min_trace(lm)) %>% 
#'   forecast()
#' 
#' @export
reconcile <- function(.data, ...){
  UseMethod("reconcile")
}

#' @rdname reconcile
#' @export
reconcile.mdl_df <- function(.data, ...){
  mutate(.data, ...)
}

#' Minimum trace forecast reconciliation
#' 
#' Reconciles a hierarchy using the minimum trace combination method. The 
#' response variable of the hierarchy must be aggregated using sums. The 
#' forecasted time points must match for all series in the hierarchy (caution:
#' this is not yet tested for beyond the series length).
#' 
#' @param models A column of models in a mable.
#' @param method The reconciliation method to use. `"mint_ridge"` is the
#'   MinT-Ridge estimator, which replaces the weight
#'   matrix `W` by `W + lambda*I` and selects the penalty `lambda` by
#'   `k`-fold cross-validation on the in-sample one-step forecasts. A penalty
#'   of `0` gives `"mint_cov"` and an unbounded penalty gives `"ols"`.
#'   Bottom-up is *not* an endpoint of this path — it corresponds to a singular
#'   reweighting that adding `lambda*I` moves away from — so use
#'   [`bottom_up()`] for that. The remaining arguments are ignored unless
#'   `method = "mint_ridge"`.
#' @param sparse If TRUE, the reconciliation will be computed using sparse 
#' matrix algebra? By default, sparse matrices will be used if the MatrixM 
#' package is installed.
#' @param k Number of cross-validation folds used to select the MinT-Ridge
#'   penalty. Must be at least 2.
#' @param window Whether the cross-validation training window expands
#'   (`"expanding"`, using all observations from the start of the series) or
#'   slides (`"rolling"`, using only the `initial` most recent observations
#'   preceding each validation block).
#' @param initial Number of in-sample observations reserved as the minimum
#'   cross-validation training window. The remaining observations are split
#'   into `k` contiguous validation blocks. If `NULL` (the default) this is
#'   `ceiling(T/2)`, capped so that at least `k` validation observations remain.
#' @param grid Penalty grid for MinT-Ridge, in units of
#'   `mean(diag(W))`. If `NULL` (the default) a `glmnet`-style data-driven grid
#'   is used: `lambda_max` is the largest eigenvalue of `W/mean(diag(W))`, and
#'   100 values are spaced geometrically from `lambda_max` down to
#'   `1e-4 * lambda_max` (`1e-2 * lambda_max` when the number of observations
#'   is not greater than the number of series), together with an exact `0` so
#'   that the unpenalised estimator remains selectable.
#' 
#' @seealso 
#' [`reconcile()`], [`aggregate_key()`]
#' 
#' @references 
#' Wickramasuriya, S. L., Athanasopoulos, G., & Hyndman, R. J. (2019). Optimal forecast reconciliation for hierarchical and grouped time series through trace minimization. Journal of the American Statistical Association, 1-45. https://doi.org/10.1080/01621459.2018.1448825 
#' 
#' @export
min_trace <- function(models, method = c("wls_var", "ols", "wls_struct", "mint_cov", "mint_shrink", "mint_ridge"),
                 sparse = NULL,
                 k = 5L,
                 window = c("expanding", "rolling"),
                 initial = NULL,
                 grid = NULL){
  if(is.null(sparse)){
    sparse <- requireNamespace("Matrix", quietly = TRUE)
  }
  method <- match.arg(method)
  window <- match.arg(window)
  if(length(k) != 1 || !is.numeric(k) || is.na(k) || k < 2 || k != round(k)){
    cli::cli_abort(c("{.arg k} must be a single number of at least 2.", "i" = "{.arg k} is {k}."))
  }
  k <- as.integer(k)
  if(!is.null(grid)){
    if(!is.numeric(grid) || length(grid) == 0 || anyNA(grid) || any(grid < 0)){
      cli::cli_abort("{.arg grid} must be a non-empty vector of non-negative, non-missing numbers.")
    }
    grid <- sort(unique(as.numeric(grid)))
  }
  if(!("mint_ridge" %in% method) &&
     (k != 5L || window != "expanding" || !is.null(initial) || !is.null(grid))){
    cli::cli_warn(c(
      "Cross-validation arguments are only used by {.code method = \"mint_ridge\"}.",
      "i" = "They are being ignored for {.code method = \"{method}\"}."
    ))
  }
  structure(models, class = c("lst_mint_mdl", "mdl_lst", "list"),
            method = method, sparse = sparse, k = k, window = window,
            initial = initial, grid = grid)
}

#' @export
forecast.lst_mint_mdl <- function(object, key_data, 
                                  new_data = NULL, h = NULL,
                                  point_forecast = list(.mean = mean), ...){
  method <- object%@%"method"
  sparse <- object%@%"sparse"
  if(sparse){
    check_installed("Matrix")
    as.matrix <- Matrix::as.matrix
    t <- Matrix::t
    diag <- function(x) if(is.vector(x)) Matrix::Diagonal(x = x) else Matrix::diag(x)
    solve <- Matrix::solve
    cov2cor <- Matrix::cov2cor
  } else {
    cov2cor <- stats::cov2cor
  }
  
  point_method <- point_forecast
  point_forecast <- list()
  # Get forecasts
  fc <- NextMethod()
  if(length(unique(map(fc, interval))) > 1){
    abort("Reconciliation of temporal hierarchies is not yet supported.")
  }
  
  # Compute weights (sample covariance)
  res <- stack_series(map(object, function(x, ...) residuals(x, ...), type = "response"))
  
  # Construct S matrix - ??GA: have moved this here as I need it for Structural scaling
  agg_data <- build_key_data_smat(key_data)
  
  n <- nrow(res)
  covm <- crossprod(stats::na.omit(res)) / n
  if(method == "ols"){
    # OLS
    W <- diag(rep(1L, nrow(covm)))
  } else if(method == "wls_var"){
    # WLS variance scaling
    W <- diag(diag(covm))
  } else if (method == "wls_struct"){
    # WLS structural scaling
    W <- diag(vapply(agg_data$agg,length,integer(1L)))
  } else if (method == "mint_cov"){
    # min_trace covariance
    W <- covm
  } else if (method == "mint_shrink"){
    # min_trace shrink
    tar <- diag(apply(res, 2, compose(crossprod, stats::na.omit))/n)
    corm <- cov2cor(covm)
    xs <- scale(res, center = FALSE, scale = sqrt(diag(covm)))
    xs <- xs[stats::complete.cases(xs),]
    v <- (1/(n * (n - 1))) * (crossprod(xs^2) - 1/n * (crossprod(xs))^2)
    diag(v) <- 0
    corapn <- cov2cor(tar)
    d <- (corm - corapn)^2
    lambda <- sum(v)/sum(d)
    lambda <- max(min(lambda, 1), 0)
    W <- lambda * tar + (1 - lambda) * covm
  } else if (method == "mint_ridge"){
    # MinT-Ridge: W + alpha*I with alpha chosen by k-fold CV on h=1 forecasts
    W <- mint_ridge_weights(covm, object, agg_data,
                            k = object%@%"k",
                            window = object%@%"window",
                            initial = object%@%"initial",
                            grid = object%@%"grid")
  } else {
    abort("Unknown reconciliation method")
  }
  
  # Check positive definiteness of weights
  if (method != "mint_ridge"){
    # mint_ridge guarantees W + alpha*I is PD by construction (alpha >= 0)
    eigenvalues <- eigen(W, only.values = TRUE)[["values"]]
    if (any(eigenvalues < 1e-8)) {
      abort("min_trace needs covariance matrix to be positive definite.", call. = FALSE)
    }
  }
  
  # Reconciliation matrices
  if(sparse){ 
    row_btm <- agg_data$leaf
    row_agg <- seq_len(nrow(key_data))[-row_btm]
    S <- Matrix::sparseMatrix(
      i = rep(seq_along(agg_data$agg), lengths(agg_data$agg)),
      j = vec_c(!!!agg_data$agg),
      x = rep(1, sum(lengths(agg_data$agg))))
    J <- Matrix::sparseMatrix(i = S[row_btm,,drop = FALSE]@i+1, j = row_btm, x = 1L, 
                              dims = rev(dim(S)))
    if (length(row_agg) == 0) {
      # Simple case of no constraints, to avoid unnecessary matrix algebra
      P <- J
    } else {
      U <- cbind(
        Matrix::Diagonal(diff(dim(J))),
        -S[row_agg,,drop = FALSE]
      )
      U <- U[, order(c(row_agg, row_btm)), drop = FALSE]
      Ut <- t(U)
      WUt <- W %*% Ut
      P <- J - J %*% WUt %*% solve(U %*% WUt, U)
      # P <- J - J%*%W%*%t(U)%*%solve(U%*%W%*%t(U))%*%U
    }
  }
  else {
    S <- build_smat_dense(agg_data)
    R <- t(S)%*%solve(W)
    P <- solve(R%*%S)%*%R
  }
  
  reconcile_fbl_list(fc, S, P, W, point_forecast = point_method)
}

# Dense summation matrix from the aggregation structure
build_smat_dense <- function(agg_data){
  S <- matrix(0L, nrow = length(agg_data$agg), ncol = max(vec_c(!!!agg_data$agg)))
  S[length(agg_data$agg)*(vec_c(!!!agg_data$agg)-1) + rep(seq_along(agg_data$agg), lengths(agg_data$agg))] <- 1L
  S
}

# Stack a list of per-series tsibbles into a time x series matrix
# Used for the residual, response and fitted values of a reconciliation model
stack_series <- function(x){
  if(length(unique(map_dbl(x, nrow))) > 1){
    # Join by index #199
    unname(as.matrix(reduce(x, full_join, by = index_var(x[[1]]))[,-1]))
  } else {
    matrix(invoke(c, map(x, `[[`, 2)), ncol = length(x))
  }
}

# glmnet-style data-driven penalty grid for MinT-Ridge.
# lambda_max is the largest eigenvalue of the scale-free W (= W / mean(diag(W))),
# the point at which alpha*I becomes commensurate with the data term and beyond
# which (W + alpha*I)^-1 is indistinguishable from alpha^-1*I, i.e. the estimator
# has converged to OLS. Note this is *not* bottom-up: the ridge path runs from
# MinT with the sample covariance (alpha = 0) to OLS (alpha = inf), and OLS is a
# different operator from bottom-up (max|P_ols - P_bu| = 1/3 for a 2-level
# hierarchy). Neither endpoint is bottom-up, which is not reachable from this
# family at all: it is the W -> diag(0,...,0,inf,...,inf) limit.
mint_ridge_grid <- function(W, Tn, n_series){
  lambda_max <- max(eigen(W, symmetric = TRUE, only.values = TRUE)[["values"]])
  # glmnet uses a wider floor when the problem is over-parameterised (T <= n)
  lambda_min <- lambda_max * if (Tn > n_series) 1e-4 else 1e-2
  c(0, exp(seq(log(lambda_min), log(lambda_max), length.out = 100)))
}

# Build the k cross-validation folds. `initial` is the minimum training window:
# the first `initial` observations are always training, and the remaining
# T - initial observations are split into k contiguous validation blocks.
# Fold j trains on rows ending at `cuts[j]` and validates on `cuts[j] + 1:cuts[j+1]`,
# so training and validation never overlap.
mint_ridge_folds <- function(Tn, k, window, initial){
  remaining <- Tn - initial
  cuts <- initial + as.integer(round(seq(0, remaining, length.out = k + 1L)))
  lapply(seq_len(k), function(j){
    list(
      train = if (window == "expanding") seq_len(cuts[j])
              else seq.int(max(cuts[j] + 1L - initial, 1L), cuts[j]),
      valid = seq.int(cuts[j] + 1L, cuts[j + 1L])
    )
  })
}

# MinT-Ridge weight matrix: W/s + alpha*I, with alpha selected by k-fold CV.
# The loss is the validation-period reconciled MSE over all series, using
# in-sample one-step (h=1) fitted values so that the covariance estimator and
# the tuning criterion are on the same footing.
mint_ridge_weights <- function(covm, object, agg_data, k, window, initial, grid){
  n_series <- nrow(covm)
  S <- build_smat_dense(agg_data)
  res <- stack_series(map(object, function(x, ...) residuals(x, ...), type = "response"))
  y   <- stack_series(map(object, response))
  fit <- stack_series(map(object, fitted))
  keep <- stats::complete.cases(res) & stats::complete.cases(y) & stats::complete.cases(fit)
  res <- res[keep, , drop = FALSE]
  y <- y[keep, , drop = FALSE]
  fit <- fit[keep, , drop = FALSE]
  Tn <- nrow(res)
  if (Tn < k + 1L){
    cli::cli_abort(c(
      "Not enough complete in-sample observations for MinT-Ridge cross-validation.",
      "i" = "{Tn} complete observation{?s} are available, but at least {k + 1L} are needed for {k} folds."
    ))
  }
  
  # Scale-free parameterisation A = W/s + alpha*I. Required: the direct
  # W + lambda*I is computationally singular for large lambda, and the scale
  # factor makes alpha invariant to the units of W.
  s <- mean(diag(covm))
  if (!is.finite(s) || s <= 0){
    cli::cli_abort(c(
      "MinT-Ridge needs a positive mean residual variance to scale the penalty.",
      "i" = "The mean diagonal of the sample covariance matrix is {s}."
    ))
  }
  W <- covm / s
  
  if (is.null(grid)){
    grid <- mint_ridge_grid(W, Tn, n_series)
  }
  if (is.null(initial)){
    initial <- min(max(ceiling(Tn / 2), 1L), Tn - k)
  }
  if (length(initial) != 1 || !is.numeric(initial) || is.na(initial) ||
      initial < 1 || initial + k > Tn){
    cli::cli_abort(c(
      "{.arg initial} must leave room for {k} cross-validation folds.",
      "i" = "There are {Tn} complete observations, so {.arg initial} must be between 1 and {Tn - k}."
    ))
  }
  initial <- as.integer(initial)
  folds <- mint_ridge_folds(Tn, k, window, initial)
  
  # alpha = 0 reproduces mint_cov, which is only defined for a positive definite
  # W. Drop it (with a warning) so that CV can pick a stabilising alpha instead.
  drop_zero <- min(eigen(W, symmetric = TRUE, only.values = TRUE)[["values"]]) <= 1e-8
  if (drop_zero && any(grid == 0)){
    cli::cli_warn(c(
      "The sample covariance matrix is not positive definite.",
      "i" = "Dropping {.code alpha = 0} from the penalty grid; cross-validation will select a stabilising penalty."
    ))
    grid <- grid[grid > 0]
  }
  
  # Per fold: eigendecompose the training covariance once, then each alpha costs
  # a single m x m solve via (W + alpha*I)^-1 = V diag(1/(d + alpha)) V'
  loss <- rep(0, length(grid))
  n_folds <- 0L
  for (fold in folds){
    tr <- fold$train
    W_tr <- crossprod(res[tr, , drop = FALSE]) / length(tr)
    s_tr <- mean(diag(W_tr))
    if (!is.finite(s_tr) || s_tr <= 0) next
    n_folds <- n_folds + 1L
    eig <- eigen(W_tr, symmetric = TRUE)
    d <- eig$values
    Q <- t(eig$vectors) %*% S
    y_va <- t(y[fold$valid, , drop = FALSE])
    f_va <- t(fit[fold$valid, , drop = FALSE])
    n_va <- length(fold$valid)
    for (i in seq_along(grid)){
      alpha <- grid[i]
      denom <- d / s_tr + alpha
      if (any(!is.finite(denom)) || any(denom <= 0)){
        loss[i] <- loss[i] + Inf
        next
      }
      D <- diag(1 / denom, nrow = length(d))
      QtD <- t(Q) %*% D
      B <- QtD %*% Q
      P <- tryCatch(solve(B) %*% QtD, error = function(e) NULL)
      if (is.null(P)){
        loss[i] <- loss[i] + Inf
      } else {
        loss[i] <- loss[i] + sum((y_va - S %*% P %*% f_va)^2) / n_va
      }
    }
  }
  
  if (n_folds == 0L){
    cli::cli_abort("MinT-Ridge cross-validation failed: no fold had a usable residual covariance.")
  }
  finite <- is.finite(loss)
  if (!any(finite)){
    cli::cli_abort("MinT-Ridge cross-validation failed: no penalty in the grid produced a usable fit.")
  }
  best <- which.min(replace(loss, !finite, Inf))
  alpha <- grid[best]
  cli::cli_inform(c(
    "i" = "MinT-Ridge: alpha = {.val {alpha}} selected by {k}-fold cross-validation ({window} window)."
  ))
  
  W + alpha * diag(n_series)
}

#' Bottom up forecast reconciliation
#' 
#' \lifecycle{experimental}
#' 
#' Reconciles a hierarchy using the bottom up reconciliation method. The 
#' response variable of the hierarchy must be aggregated using sums. The 
#' forecasted time points must match for all series in the hierarchy.
#' 
#' @param models A column of models in a mable.
#' 
#' @seealso 
#' [`reconcile()`], [`aggregate_key()`]
#' @export
bottom_up <- function(models){
  structure(models, class = c("lst_btmup_mdl", "mdl_lst", "list"))
}

#' @export
forecast.lst_btmup_mdl <- function(object, key_data, 
                                   point_forecast = list(.mean = mean),
                                   new_data = NULL, ...){
  # Keep only bottom layer
  agg_data <- build_key_data_smat(key_data)
  
  S <- matrix(0L, nrow = length(agg_data$agg), ncol = max(vec_c(!!!agg_data$agg)))
  S[length(agg_data$agg)*(vec_c(!!!agg_data$agg)-1) + rep(seq_along(agg_data$agg), lengths(agg_data$agg))] <- 1L
  
  btm <- agg_data$leaf
  # object <- object[btm]
  # if(!is.null(new_data)){
  #   new_data <- new_data[btm]
  # }
  
  point_method <- point_forecast
  point_forecast <- list()
  
  # Get base forecasts
  # fc <- vector("list", nrow(S))
  # fc[btm] <- NextMethod()
  fc <- NextMethod()
  
  # Add dummy forecasts to unused levels
  # fc[seq_along(fc)[-btm]] <- fc[btm[1]]
  
  P <- matrix(0L, nrow = ncol(S), ncol = nrow(S))
  P[(btm-1L)*nrow(P) + seq_len(nrow(P))] <- 1L
  
  reconcile_fbl_list(fc, S, P, W = diag(nrow(S)),
                     point_forecast = point_method)
}


#' Top down forecast reconciliation
#' 
#' \lifecycle{experimental}
#' 
#' Reconciles a hierarchy using the top down reconciliation method. The 
#' response variable of the hierarchy must be aggregated using sums. The 
#' forecasted time points must match for all series in the hierarchy.
#' 
#' @param models A column of models in a mable.
#' @param method The reconciliation method to use.
#' 
#' @seealso 
#' [`reconcile()`], [`aggregate_key()`]
#' 
#' @export
top_down <- function(models, method = c("forecast_proportions", "average_proportions", "proportion_averages")){
  structure(models, class = c("lst_topdwn_mdl", "mdl_lst", "list"),
            method = match.arg(method))
}

#' @export
forecast.lst_topdwn_mdl <- function(object, key_data, 
                                    point_forecast = list(.mean = mean), ...){
  method <- object%@%"method"
  point_method <- point_forecast
  point_forecast <- list()
  
  agg_data <- build_key_data_smat(key_data)
  S <- matrix(0L, nrow = length(agg_data$agg), ncol = max(vec_c(!!!agg_data$agg)))
  S[length(agg_data$agg)*(vec_c(!!!agg_data$agg)-1) + rep(seq_along(agg_data$agg), lengths(agg_data$agg))] <- 1L
  # Identify top and bottom level
  top <- which.max(rowSums(S))
  btm <- agg_data$leaf
  
  kv <- names(key_data)[-ncol(key_data)]
  agg_shadow <- as_tibble(map(key_data[kv], is_aggregated))
  agg_struct <- vctrs::vec_unique(agg_shadow)
  agg_depth <- nrow(agg_struct)
  if(length(kv) != (agg_depth - 1)) {
    abort("Top down reconciliation requires strictly hierarchical structures.")
  }
  agg_order <- kv[order(vapply(agg_struct, sum, integer(1L)))]
  
  if(method == "forecast_proportions") {
    fc <- NextMethod()
    fc_dist <- lapply(fc, function(x) x[[distribution_var(x)]])
    fc_mean <- lapply(fc_dist, mean)
    fc_mean <- do.call(cbind, fc_mean)
    # Ensure key structure matches order of fc object
    key_data <- key_data[order(vec_c(!!!key_data$.rows)),]
    level_id <- rowSums(agg_shadow)
    td <- propagate_forecast_proportions(
      fc_mean, key_data, agg_order, level_id,
      start_nodes = which(level_id == length(kv))
    )
    fc_prop <- td$fc_prop
    # Code adapted from reconcile_fbl_list to handle changing weights over horizon
    # This will need to be refactored later so that reconcile_fbl_list is broken up into more sub-problems
    # As the weight matrix is an identity, this code and computation is much simpler.
    is_normal <- all(map_lgl(fc_dist, function(x) all(dist_types(x) == "dist_normal")))
    # Point forecast means can be computed in one step
    fc_mean <- split(fc_mean[,top]*fc_prop,col(fc_prop))
    if(is_normal) {
      fc_var <- map(fc_dist, distributional::variance)
      fc_var <- fc_prop * fc_var[[top]] * fc_prop
      fc_var <- split(fc_var, col(fc_var))
      fc_dist <- map2(fc_mean, map(fc_var, sqrt), distributional::dist_normal)
    } else {
      fc_dist <- lapply(fc_mean, distributional::dist_degenerate)
    }
    # Update fables
    fc <- map2(fc, fc_dist, function(fc, dist){
      dimnames(dist) <- dimnames(fc[[distribution_var(fc)]])
      fc[[distribution_var(fc)]] <- dist
      point_fc <- compute_point_forecasts(dist, point_method)
      fc[names(point_fc)] <- point_fc
      fc
    })
    return(fc)
    
  } else {
    # Compute dis-aggregation matrix
    history <- lapply(object, function(x) response(x)[[".response"]])
    top_y <- history[[top]]
    btm_y <- history[btm]
    if (method == "average_proportions") { 
      prop <- map_dbl(btm_y, function(y) mean(y/top_y))
    } else if (method == "proportion_averages") {
      prop <- map_dbl(btm_y, mean) / mean(top_y)
    } else {
      abort("Unkown `top_down()` reconciliation `method`.")
    }
    
    # Keep only top layer
    object <- object[top]
    
    # Get base forecasts
    fc <- vector("list", nrow(S))
    fc[top] <- NextMethod()
    
    # Add dummy forecasts to unused levels
    fc[seq_along(fc)[-top]] <- fc[top]
  }
  
  P <- matrix(0L, nrow = ncol(S), ncol = nrow(S))
  P[,top] <- prop
  
  reconcile_fbl_list(fc, S, P, W = diag(nrow(S)),
                     point_forecast = point_method)
}


# Propagate top-down forecast proportions from a set of start nodes to all
# descendants, one level at a time.
#
# Returns a list with:
#   fc_prop  - (h x n_all) matrix: 1 at start_nodes, P(j | root_j) for
#              descendants, 0 for nodes above the start level
#   root_idx - integer vector (n_all): position within start_nodes for each
#              node's root; 0 for nodes not descended from start_nodes
propagate_forecast_proportions <- function(fc_mean, key_data, agg_order,
                                           level_id, start_nodes) {
  h        <- nrow(fc_mean)
  n_all    <- ncol(fc_mean)
  n_start  <- length(start_nodes)
  start_level <- level_id[[start_nodes[[1L]]]]

  fc_prop <- matrix(0, h, n_all)
  fc_prop[, start_nodes] <- 1
  root_idx <- integer(n_all)
  root_idx[start_nodes] <- seq_len(n_start)

  # Compact (h x current-layer) proportions, replaced at each level
  fc_prop_layer <- matrix(1, h, n_start)
  parent_loc <- start_nodes

  for (i in seq.int(length(agg_order) - start_level + 1L, length(agg_order))) {
    child_loc  <- which(level_id == (length(agg_order) - i))
    agg_parent <- vec_match(
      key_data[child_loc,  agg_order[seq_len(i - 1L)], drop = FALSE],
      key_data[parent_loc, agg_order[seq_len(i - 1L)], drop = FALSE]
    )
    if (anyNA(agg_parent)) {
      abort("An error has occurred when reconciling the hierarchy.\nPlease report this bug here: https://github.com/tidyverts/fabletools/issues")
    }
    fc_prop_layer <- fc_prop_layer[, agg_parent, drop = FALSE] *
      fc_mean[, child_loc, drop = FALSE] /
      t(rowsum(t(fc_mean[, child_loc, drop = FALSE]), agg_parent))[, agg_parent]
    fc_prop[, child_loc] <- fc_prop_layer
    root_idx[child_loc]  <- root_idx[parent_loc][agg_parent]
    parent_loc <- child_loc
  }

  list(fc_prop = fc_prop, root_idx = root_idx)
}

#' Middle out forecast reconciliation
#' 
#' \lifecycle{experimental}
#' 
#' Reconciles a hierarchy using the middle out reconciliation method. The 
#' response variable of the hierarchy must be aggregated using sums. The 
#' forecasted time points must match for all series in the hierarchy.
#' 
#' @param models A column of models in a mable.
#' @param split The middle level of the hierarchy from which the bottom-up and
#' top-down approaches are used above and below respectively.
#' 
#' @seealso 
#' [`reconcile()`], [`aggregate_key()`]
#' [*Forecasting: Principles and Practice* - Middle-out approach](https://otexts.com/fpp3/single-level.html#middle-out-approach)
#' 
#' @export
middle_out <- function(models, split = 1){
  structure(models, class = c("lst_midout_mdl", "mdl_lst", "list"),
            split = split)
}

#' @export
forecast.lst_midout_mdl <- function(object, key_data, 
                                    point_forecast = list(.mean = mean),
                                    new_data = NULL, ...){
  split <- object%@%"split"
  point_method <- point_forecast
  point_forecast <- list()
  
  agg_data <- build_key_data_smat(key_data)
  S <- matrix(0L, nrow = length(agg_data$agg), ncol = max(vec_c(!!!agg_data$agg)))
  S[length(agg_data$agg)*(vec_c(!!!agg_data$agg)-1) + rep(seq_along(agg_data$agg), lengths(agg_data$agg))] <- 1L
  
  # Identify top and bottom level
  top <- which.max(rowSums(S))
  btm <- agg_data$leaf
  
  key_data <- key_data[order(vec_c(!!!key_data$.rows)),]
  object <- object[order(vec_c(!!!key_data$.rows))]
  if(!is.null(new_data)){
    new_data <- new_data[order(vec_c(!!!key_data$.rows)),]
  }
  
  kv <- names(key_data)[-ncol(key_data)]
  agg_shadow <- as_tibble(map(key_data[kv], is_aggregated))
  agg_struct <- vctrs::vec_unique(agg_shadow)
  agg_depth <- nrow(agg_struct)
  if(length(kv) != (agg_depth - 1)) {
    abort("Middle out reconciliation requires strictly hierarchical structures.")
  }
  agg_order <- kv[order(vapply(agg_struct, sum, integer(1L)))]
  if(is.character(split)) {
    split <- match(split, agg_order)
  }
  if(is.na(split) || split < 1L || split > length(agg_order)) {
    abort("`split` must identify one of the hierarchy levels.")
  }
  
  level_id <- rowSums(agg_shadow)
  split_level <- length(agg_order) - split
  split_nodes <- which(level_id == split_level)
  nodes_above <- which(level_id > split_level)
  object <- object[-nodes_above]
  if(!is.null(new_data)){
    new_data <- new_data[-nodes_above]
  }
  fc <- NextMethod()
  
  fc_dist <- lapply(fc, function(x) x[[distribution_var(x)]])
  h <- vec_size(fc_dist[[1]])
  fc_mean <- matrix(0, h, nrow(key_data))
  fc_mean[,-nodes_above] <- do.call(cbind, lapply(fc_dist, mean))
  
  td <- propagate_forecast_proportions(
    fc_mean, key_data, agg_order, level_id,
    start_nodes = split_nodes
  )
  btm_loc        <- which(level_id == 0L)
  mid_root_nodes <- split_nodes[td$root_idx[btm_loc]]
  fc_prop        <- td$fc_prop[, btm_loc, drop = FALSE]

  fc_mean <- (fc_prop * fc_mean[, mid_root_nodes]) %*% t(S)
  # Code adapted from reconcile_fbl_list to handle changing weights over horizon
  # This will need to be refactored later so that reconcile_fbl_list is broken up into more sub-problems
  # As the weight matrix is an identity, this code and computation is much simpler.
  is_normal <- all(map_lgl(fc_dist, function(x) all(dist_types(x) == "dist_normal")))
  # Point forecast means can be computed in one step
  fc_mean <- split(fc_mean,col(fc_mean))
  if(is_normal) {
    fc_var <- vector("list", nrow(key_data))
    fc_var[-nodes_above] <- map(fc_dist, distributional::variance)
    fc_var[nodes_above] <- rep_len(list(double(h)), length(nodes_above))
    
    P <- matrix(0L, nrow = ncol(S), ncol = nrow(S))
    
    # (S%*%P)%*%t(fc_mean)
    fc_var <- map(seq_len(h), function(i) {
      # Add top down structure
      P[seq_along(mid_root_nodes) + (mid_root_nodes-1)*nrow(P)] <- fc_prop[i,]
      SP <- S%*%P
      diag(SP%*%diag(map_dbl(fc_var, `[[`, i))%*%t(SP))
    })
    fc_dist <- map2(fc_mean, transpose_dbl(map(fc_var, sqrt)), distributional::dist_normal)
    
  } else {
    fc_dist <- lapply(fc_mean, distributional::dist_degenerate)
  }
  
  # Update fables
  map2(rep(fc[1], nrow(key_data)), fc_dist, function(fc, dist){
    dimnames(dist) <- dimnames(fc[[distribution_var(fc)]])
    fc[[distribution_var(fc)]] <- dist
    point_fc <- compute_point_forecasts(dist, point_method)
    fc[names(point_fc)] <- point_fc
    fc
  })
}

reconcile_fbl_list <- function(fc, S, P, W, point_forecast, SP = NULL) {
  if(length(unique(map(fc, interval))) > 1){
    abort("Reconciliation of temporal hierarchies is not yet supported.")
  }
  if(!inherits(S, "matrix")) {
    # Use sparse functions
    check_installed("Matrix")
    as.matrix <- Matrix::as.matrix
    t <- Matrix::t
    diag <- function(x) if(is.vector(x)) Matrix::Diagonal(x = x) else Matrix::diag(x)
    cov2cor <- Matrix::cov2cor
  } else {
    cov2cor <- stats::cov2cor
  }
  if(is.null(SP)) {
    SP <- S%*%P
  }
  
  fc_dist <- map(fc, function(x) x[[distribution_var(x)]])
  dist_type <- lapply(fc_dist, function(x) unique(dist_types(x)))
  dist_type <- unique(unlist(dist_type))
  is_normal <- all(map_lgl(fc_dist, function(x) all(dist_types(x) == "dist_normal")))
  
  fc_mean <- as.matrix(invoke(cbind, map(fc_dist, mean)))
  fc_var <- transpose_dbl(map(fc_dist, distributional::variance))
  
  # Apply to forecasts
  fc_mean <- as.matrix(SP%*%t(fc_mean))
  fc_mean <- split(fc_mean, row(fc_mean))
  if(identical(dist_type, "dist_normal")){
    R1 <- cov2cor(W)
    W_h <- map(fc_var, function(var) diag(sqrt(var))%*%R1%*%t(diag(sqrt(var))))
    fc_var <- map(W_h, function(W) diag(SP%*%W%*%t(SP)))
    fc_dist <- map2(fc_mean, transpose_dbl(map(fc_var, sqrt)), distributional::dist_normal)
  } else if (identical(dist_type, "dist_sample")) {
    sample_size <- unique(unlist(lapply(fc_dist, function(x) unique(lengths(distributional::parameters(x)$x)))))
    if(length(sample_size) != 1L) stop("Cannot reconcile sample paths with different replication sizes.")
    sample_horizon <- unique(lengths(fc_dist))
    if(length(sample_horizon) != 1L) stop("Cannot reconcile sample paths with different forecast horizon lengths.")
    # Extract sample paths
    samples <- lapply(fc_dist, function(x) distributional::parameters(x)$x)
    # Convert to array [samples,horizon,nodes]
    samples <- array(unlist(samples, use.names = FALSE), dim = c(sample_size, sample_horizon, length(fc_dist)))
    # Reconcile
    samples <- apply(samples, 1, function(x) as.matrix(SP%*%t(x)), simplify = FALSE)
    # Convert to array [nodes, horizon, samples]
    samples <- array(unlist(samples), dim = c(length(fc_dist), sample_horizon, sample_size))
    # Convert to distributions
    fc_dist <- apply(
      samples, 1L, simplify = FALSE,
      function(x) unname(distributional::dist_sample(split(x, row(x))))
    )
  } else {
    fc_dist <- map(fc_mean, distributional::dist_degenerate)
  }
  
  # Update fables
  map2(fc, fc_dist, function(fc, dist){
    dimnames(dist) <- dimnames(fc[[distribution_var(fc)]])
    fc[[distribution_var(fc)]] <- dist
    point_fc <- compute_point_forecasts(dist, point_forecast)
    fc[names(point_fc)] <- point_fc
    fc
  })
}

build_smat_rows <- function(key_data){
  lifecycle::deprecate_warn("0.2.1", "fabletools::build_smat_rows()", "fabletools::build_key_data_smat()")
  row_col <- sym(colnames(key_data)[length(key_data)])
  
  smat <- key_data %>%
    unnest(!!row_col) %>% 
    dplyr::arrange(!!row_col) %>% 
    select(!!expr(-!!row_col))
  
  agg_struc <- group_data(dplyr::group_by_all(as_tibble(map(smat, is_aggregated))))
  
  # key_unique <- map(smat, function(x){
  #   x <- unique(x)
  #   x[!is_aggregated(x)]
  # })
  
  agg_struc$.smat <- map(agg_struc$.rows, function(n) diag(1, nrow = length(n), ncol = length(n)))
  agg_struc <- map(seq_len(nrow(agg_struc)), function(i) agg_struc[i,])
  
  out <- reduce(agg_struc, function(x, y){
    # For now, assume x is aggregated into y somehow
    n_key <- ncol(x)-2
    nm_key <- names(x)[seq_len(n_key)]
    agg_vars <- map2_lgl(x[seq_len(n_key)], y[seq_len(n_key)], `<`)
    
    if(!any(agg_vars)) abort("Something unexpected happened, please report this bug at https://github.com/tidyverts/fabletools/issues/ with a description of what you're trying to do.")
    
    # Match rows between summation matrices
    not_agg <- names(Filter(`!`, y[seq_len(n_key)]))
    cols <- group_data(group_by(smat[x$.rows[[1]][seq_len(ncol(x$.smat[[1]]))],], !!!syms(not_agg)))$.rows
    cols_pos <- unlist(cols)
    cols <- rep(seq_along(cols), map_dbl(cols, length))
    cols[cols_pos] <- cols
    
    x$.rows[[1]] <- c(x$.rows[[1]], y$.rows[[1]])
    x$.smat <- list(rbind(
      x$.smat[[1]],
      y$.smat[[1]][, cols, drop = FALSE]
    ))
    x
  })
  
  smat <- out$.smat[[1]]
  smat[out$.rows[[1]],] <- smat
  
  return(smat)
}

build_key_data_smat <- function(x){
  if (any(lengths(x[[ncol(x)]]) > 1L)) {
    # Find the order based on the first entry position of each key
    x[[ncol(x)]] <- as.list(rank(vapply(x[[ncol(x)]], min, integer(1L))))
  }

  kv <- names(x)[-ncol(x)]
  agg_shadow <- as_tibble(map(x[kv], is_aggregated))
  grp <- as_tibble(vctrs::vec_group_loc(agg_shadow))
  num_agg <- rowSums(grp$key)
  # Initialise comparison leafs with known/guaranteed leafs
  x_leaf <- x[unlist(grp$loc[which(num_agg == min(num_agg))]),]
  
  # Sort by disaggregation to identify aggregated leafs in order
  grp <- grp[order(num_agg),]
  
  grp$match <- lapply(unname(split(grp, seq_len(nrow(grp)))), function(level){
    disagg_col <- which(!vec_c(!!!level$key))
    agg_idx <- level[["loc"]][[1]]
    pos <- vec_match(x_leaf[disagg_col], x[agg_idx, disagg_col])
    pos <- vec_group_loc(pos)
    pos <- pos[!is.na(pos$key),]
    # Add non-matches as leaf nodes
    agg_leaf <- setdiff(seq_along(agg_idx), pos$key)
    if(!is_empty(agg_leaf)){
      pos <- vec_rbind(
        pos,
        structure(list(key = agg_leaf, loc = as.list(seq_along(agg_leaf) + nrow(x_leaf))), 
                  class = "data.frame", row.names = agg_leaf)
      )
      x_leaf <<- vec_rbind(
        x_leaf, 
        x[agg_idx[agg_leaf],]
      )
    }
    pos$loc[order(pos$key)]
  })
  if(any(lengths(grp$loc) != lengths(grp$match))) {
    abort("An error has occurred when constructing the summation matrix.\nPlease report this bug here: https://github.com/tidyverts/fabletools/issues")
  }
  idx_leaf <- unlist(x_leaf$.rows)
  x$.rows[unlist(x$.rows)[unlist(grp$loc)]] <- unlist(grp$match, recursive = FALSE)
  return(list(agg = x$.rows, leaf = idx_leaf))
}
