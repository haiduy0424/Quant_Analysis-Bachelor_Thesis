# A. =========================================================================== ------------ ENGINEER ------------ =========================================================================  

# 0. SETUP ==================================================================================================================================================================================

library(roll)
library(xts)
library(dccmidas)
library(rumidas)
library(readxl)
library(quantmod)
library(dplyr)
library(openxlsx)
library(tseries)
library(strucchange)
library(ggplot2)
library(lubridate)
library(forecast)
library(maxLik)     
library(numDeriv)   

# 1. UTILITIES & STRUCTURAL BREAK TEST FUNCTION =============================================================================================================================================

# Create wrapper for dccmidas:::QMLE_sd 
QMLE_sd_pkg <- dccmidas:::QMLE_sd

qmle_inference_engine <- function(est, ll_func = NULL, ll_args = NULL,
                         daily_ret = NULL, mv_m = NULL, mv_m_m = NULL,
                         K = NULL, dummy = NULL, n_boot = 200) {

  n_params <- length(coef(est))
  param_names <- names(coef(est))
  
  se_qmle <- tryCatch(QMLE_sd_pkg(est), error = function(e) NULL)
  if (!is.null(se_qmle) && all(!is.na(se_qmle)) && all(is.finite(se_qmle))) {
    cat("  [SE] QMLE sandwich SE: OK\n")
    return(se_qmle)
  }
  cat("  [SE] QMLE sandwich failed.\n")

  if (!is.null(ll_func) && !is.null(ll_args)) {

    theta_hat <- coef(est)
    boundary_flags <- c()
    if ("alpha" %in% param_names && theta_hat["alpha"] < 0.005)
      boundary_flags <- c(boundary_flags, "alpha near 0")
    if ("beta" %in% param_names && theta_hat["beta"] > 0.995)
      boundary_flags <- c(boundary_flags, "beta near 1")
    if (any(grepl("^w2", param_names))) {
      w2_params <- theta_hat[grepl("^w2", param_names)]
      if (any(w2_params < 1.01))
        boundary_flags <- c(boundary_flags, paste("omega near boundary:",
                                                  paste(names(w2_params[w2_params < 1.01]), collapse=", ")))
    }
    if (length(boundary_flags) > 0) {
      cat("  [SE] WARNING — boundary solutions detected:\n")
      for (f in boundary_flags) cat("        ", f, "\n")
      cat("        SEs may be unreliable. Consider re-estimating with wider grid.\n")
    }

    sum_ll <- function(p) {
      args <- c(list(param = p), ll_args)
      sum(do.call(ll_func, args))
    }

    H <- tryCatch({
      -numDeriv::hessian(func = sum_ll, x = coef(est))
    }, error = function(e) NULL)

    if (!is.null(H) && all(is.finite(H))) {
   
      cond_num <- tryCatch(kappa(H, exact = TRUE), error = function(e) Inf)
      if (cond_num < 1e12) {
        inv_H <- tryCatch(solve(H), error = function(e) NULL)
        if (!is.null(inv_H)) {
          se_diag <- diag(inv_H)
          se_diag[se_diag < 0] <- NA
          se_hess <- sqrt(se_diag)
          if (any(!is.na(se_hess))) {
            cat("  [SE] Numerical Hessian: OK (condition number:",
                formatC(cond_num, format = "e", digits = 2), ")\n")
            return(se_hess)
          }
        }
      } else {
        cat("  [SE] Hessian near-singular (condition:",
            formatC(cond_num, format = "e", digits = 2), "). Skipping.\n")
      }
    } else {
      cat("  [SE] Numerical Hessian contains Inf/NaN. Skipping.\n")
    }
  }

  if (!is.null(ll_func) && !is.null(ll_args)) {
    cat("  [SE] Trying OPG sandwich estimator...\n")

    se_opg <- tryCatch({
      theta_hat <- coef(est)

      G <- numDeriv::jacobian(func = function(p) {
        args <- c(list(param = p), ll_args)
        do.call(ll_func, args) 
      }, x = theta_hat)

      B <- crossprod(G)  #

      sum_ll_local <- function(p) {
        args <- c(list(param = p), ll_args)
        sum(do.call(ll_func, args))
      }
      A <- tryCatch(-numDeriv::hessian(sum_ll_local, theta_hat), error = function(e) NULL)

      if (!is.null(A) && all(is.finite(A)) && kappa(A) < 1e12) {
        A_inv <- solve(A)
        V_sandwich <- A_inv %*% B %*% A_inv
        se_vals <- sqrt(pmax(diag(V_sandwich), 0))
        cat("  [SE] OPG sandwich: OK\n")
        se_vals
      } else {
        V_opg <- solve(B)
        se_vals <- sqrt(pmax(diag(V_opg), 0))
        cat("  [SE] OPG (no sandwich, Hessian unavailable): OK\n")
        se_vals
      }
    }, error = function(e) {
      cat("  [SE] OPG failed:", e$message, "\n")
      NULL
    })

    if (!is.null(se_opg) && any(!is.na(se_opg))) return(se_opg)
  }

  if (!is.null(ll_func) && !is.null(ll_args) && !is.null(daily_ret)) {
    cat("  [SE] Trying parametric bootstrap (n=", n_boot, ")...\n")

    se_boot <- tryCatch({
      theta_hat <- coef(est)
      TT <- length(daily_ret)
      boot_params <- matrix(NA, nrow = n_boot, ncol = n_params)

      orig_ll <- do.call(ll_func, c(list(param = theta_hat), ll_args))

      successful <- 0
      attempts <- 0
      max_attempts <- n_boot * 3

      while (successful < n_boot && attempts < max_attempts) {
        attempts <- attempts + 1

        block_size <- max(5L, round(TT^(1/3)))
        n_blocks <- ceiling(TT / block_size)
        block_starts <- sample(1:(TT - block_size + 1), n_blocks, replace = TRUE)
        indices <- unlist(lapply(block_starts, function(s) s:(s + block_size - 1)))
        indices <- indices[1:TT]

        boot_ret <- daily_ret[indices]

        boot_ll_args <- ll_args
        boot_ll_args$daily_ret <- as.numeric(boot_ret)

        boot_est <- tryCatch({
          suppressWarnings(maxLik(
            logLik = ll_func, start = theta_hat,
            daily_ret = as.numeric(boot_ret),
            mv_m = ll_args$mv_m,
            mv_m_m = if ("mv_m_m" %in% names(ll_args)) ll_args$mv_m_m else NULL,
            K = ll_args$K,
            distribution = if ("distribution" %in% names(ll_args)) ll_args$distribution else "norm",
            lag_fun = if ("lag_fun" %in% names(ll_args)) ll_args$lag_fun else "Beta",
            dummy = if ("dummy" %in% names(ll_args)) ll_args$dummy else NULL,
            constraints = list(
              ineqA = rbind(c(1, rep(0, n_params-1)),
                            c(0, 1, rep(0, n_params-2)),
                            c(-1, -1, rep(0, n_params-2))),
              ineqB = c(-0.001, -0.001, 0.999)
            ),
            iterlim = 500, method = "BFGS"
          ))
        }, error = function(e) NULL)

        if (!is.null(boot_est) && boot_est$code %in% c(0, 1, 2)) {
          successful <- successful + 1
          boot_params[successful, ] <- coef(boot_est)
          if (successful %% 50 == 0) cat("    Bootstrap:", successful, "/", n_boot, "\n")
        }
      }

      cat("    Bootstrap completed:", successful, "successful out of", attempts, "attempts\n")

      if (successful >= 30) {
        se_vals <- apply(boot_params[1:successful, ], 2, sd, na.rm = TRUE)
        cat("  [SE] Bootstrap SE: OK\n")
        se_vals
      } else {
        cat("  [SE] Bootstrap: too few successful replications.\n")
        NULL
      }
    }, error = function(e) {
      cat("  [SE] Bootstrap failed:", e$message, "\n")
      NULL
    })

    if (!is.null(se_boot) && any(!is.na(se_boot))) return(se_boot)
  }

  if (!is.null(ll_func) && !is.null(ll_args)) {
    cat("  [SE] Trying profile likelihood for individual parameters...\n")

    se_profile <- tryCatch({
      theta_hat <- coef(est)
      max_ll <- sum(do.call(ll_func, c(list(param = theta_hat), ll_args)))
      se_vals <- rep(NA, n_params)

      for (j in 1:n_params) {
        target <- max_ll - 1.92

        profile_ll <- function(pj) {
          opt_func <- function(p_rest) {
            p_full <- theta_hat
            p_full[j] <- pj
            p_full[-j] <- p_rest
            -sum(do.call(ll_func, c(list(param = p_full), ll_args)))
          }

          opt <- tryCatch(
            optim(theta_hat[-j], opt_func, method = "Nelder-Mead",
                  control = list(maxit = 200)),
            error = function(e) NULL
          )
          if (is.null(opt)) return(-Inf)
          return(-opt$value)
        }

        upper <- tryCatch({
          uniroot(function(pj) profile_ll(pj) - target,
                  interval = c(theta_hat[j], theta_hat[j] + 5 * abs(theta_hat[j]) + 1),
                  extendInt = "upX", maxiter = 30)$root
        }, error = function(e) NA)

        lower <- tryCatch({
          uniroot(function(pj) profile_ll(pj) - target,
                  interval = c(theta_hat[j] - 5 * abs(theta_hat[j]) - 1, theta_hat[j]),
                  extendInt = "downX", maxiter = 30)$root
        }, error = function(e) NA)

        if (!is.na(upper) && !is.na(lower)) {
          se_vals[j] <- (upper - lower) / (2 * 1.96)
          cat("    param", j, "(", param_names[j], "): profile CI [",
              round(lower, 4), ",", round(upper, 4), "] -> SE =", round(se_vals[j], 4), "\n")
        }
      }

      if (any(!is.na(se_vals))) {
        cat("  [SE] Profile likelihood: OK (", sum(!is.na(se_vals)), "/", n_params, "params)\n")
        se_vals
      } else {
        NULL
      }
    }, error = function(e) {
      cat("  [SE] Profile likelihood failed:", e$message, "\n")
      NULL
    })

    if (!is.null(se_profile) && any(!is.na(se_profile))) return(se_profile)
  }

  cat("  [SE] ALL METHODS FAILED. Returning NAs.\n")
  cat("  [SE] This typically indicates:\n")
  cat("        - Flat likelihood surface (weak identification)\n")
  cat("        - Parameter at boundary of parameter space\n")
  cat("        - Model misspecification for this asset\n")
  return(rep(NA, n_params))
}

# Bai-Perron structural break test & Dummy variable construction
bai_perron_test <- function(series, max_breaks = NULL,
                            start_date = NULL, end_date = NULL,
                            selection = "BIC", min_supF_pval = 0.05) {
  if (!inherits(series, "xts")) stop("Error: Input must be an xts time series.")

  if (!is.null(start_date) && !is.null(end_date)) {
    series <- series[paste(start_date, end_date, sep = "/")]
  }

  numeric_values <- as.numeric(series)
  time_index <- index(series)

  if (any(is.na(numeric_values)) || length(numeric_values) < 10) {
    stop("Error: The series contains NAs or is too short for analysis.")
  }
  
  if (is.null(max_breaks)) {
    bp_result <- breakpoints(numeric_values ~ 1)
    max_candidates <- length(bp_result$breakpoints)
  } else {
    bp_result <- breakpoints(numeric_values ~ 1, breaks = max_breaks)
    max_candidates <- max_breaks
  }

  cat("\n", strrep("=", 60), "\n")
  cat(" Bai-Perron Structural Break Test (data-driven)\n")
  cat(strrep("=", 60), "\n")
  cat(sprintf("Sample size T = %d, max candidate breaks = %d\n",
              length(numeric_values), max_candidates))

  bic_vals <- numeric(max_candidates + 1)
  bic_vals[1] <- BIC(lm(numeric_values ~ 1))  

  for (nb in 1:max_candidates) {
    bp_nb <- tryCatch(breakpoints(bp_result, breaks = nb),
                      error = function(e) NULL)
    if (!is.null(bp_nb)) {
      bic_vals[nb + 1] <- AIC(bp_nb, k = log(length(numeric_values)))
    } else {
      bic_vals[nb + 1] <- Inf
    }
  }

  bic_optimal <- which.min(bic_vals) - 1L   

  cat("\nBIC values by # of breaks:\n")
  for (nb in 0:max_candidates) {
    flag <- if (nb == bic_optimal) "  <-- BIC OPTIMAL" else ""
    cat(sprintf("  %d breaks: BIC = %.4f%s\n", nb, bic_vals[nb + 1], flag))
  }

  fs_test <- tryCatch(sctest(Fstats(numeric_values ~ 1), type = "supF"),
                      error = function(e) NULL)
  supF_significant <- !is.null(fs_test) && fs_test$p.value < min_supF_pval

  cat("\nsupF test (null: no breaks):\n")
  if (!is.null(fs_test)) {
    cat(sprintf("  supF statistic = %.4f, p-value = %.6f\n",
                fs_test$statistic, fs_test$p.value))
    cat(sprintf("  -> %s null at %.0f%% level\n",
                if (supF_significant) "REJECT" else "FAIL TO REJECT",
                min_supF_pval * 100))
  } else {
    cat("  supF test failed to run.\n")
  }

  if (selection == "BIC") {
    n_breaks_final <- bic_optimal
  } else if (selection == "supF") {
    n_breaks_final <- if (supF_significant) max_candidates else 0L
  } else {
    stop("selection must be 'BIC' or 'supF'")
  }

  if (n_breaks_final > 0 && !supF_significant) {
    cat(sprintf("\nNOTE: BIC suggests %d break(s), but supF cannot reject 'no breaks' at %.0f%%.\n",
                n_breaks_final, min_supF_pval * 100))
    cat("      Defaulting to zero breaks (parsimony).\n")
    n_breaks_final <- 0L
  }

  if (n_breaks_final == 0) {
    cat("\nFINAL DECISION: 0 breaks selected. SB models will reduce to no-SB models.\n")
    cat(strrep("=", 60), "\n\n")
    return(list(breakpoints_vector       = NULL,
                breakpoints_list         = list(),
                breakpoints_index_vector = integer(0),
                breakpoints_index_list   = list(),
                n_breaks                 = 0L,
                bic_optimal              = bic_optimal,
                bic_values               = bic_vals,
                supF_test                = fs_test,
                confidence_intervals     = NULL))
  }

  bp_final <- breakpoints(bp_result, breaks = n_breaks_final)
  breakpoints_index <- bp_final$breakpoints
  breakpoint_dates <- time_index[breakpoints_index]

  ci_bp <- tryCatch(confint(bp_result, breaks = n_breaks_final),
                    error = function(e) NULL)

  if (!is.null(ci_bp)) {
    cat("\nBreakpoint confidence intervals (95%):\n")
    for (j in 1:n_breaks_final) {
      lo  <- time_index[ci_bp$confint[j, 1]]
      mid <- time_index[ci_bp$confint[j, 2]]
      hi  <- time_index[ci_bp$confint[j, 3]]
      cat(sprintf("  Break %d: %s  [%s, %s]\n", j, mid, lo, hi))
    }
  }

  cat(sprintf("\nFINAL DECISION: %d break(s) selected at: %s\n",
              n_breaks_final, paste(breakpoint_dates, collapse = ", ")))
  cat(strrep("=", 60), "\n\n")

  # Plot 
  df <- data.frame(Time = time_index, Value = numeric_values)
  p <- ggplot(df, aes(x = Time, y = Value)) +
    geom_line(color = "blue", size = 1) +
    geom_vline(xintercept = as.numeric(breakpoint_dates),
               color = "red", linetype = "dashed", size = 1) +
    labs(title = sprintf("Bai-Perron Test: %d break(s) detected", n_breaks_final),
         subtitle = paste("Breakpoints at:", paste(breakpoint_dates, collapse = ", ")),
         x = "Time", y = "Series Level") +
    theme_minimal()

  if (!is.null(ci_bp)) {
    for (j in 1:n_breaks_final) {
      p <- p + annotate("rect",
                        xmin = as.numeric(time_index[ci_bp$confint[j, 1]]),
                        xmax = as.numeric(time_index[ci_bp$confint[j, 3]]),
                        ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "red")
    }
  }
  print(p)

  return(list(breakpoints_vector       = breakpoint_dates,
              breakpoints_list         = as.list(breakpoint_dates),
              breakpoints_index_vector = breakpoints_index,
              breakpoints_index_list   = as.list(breakpoints_index),
              n_breaks                 = n_breaks_final,
              bic_optimal              = bic_optimal,
              bic_values               = bic_vals,
              supF_test                = fs_test,
              confidence_intervals     = ci_bp))
}

# Create dummy variables
create_dummy_variables <- function(breakpoints_dates, series) {
  library(xts)

  if (!inherits(series, "xts")) {
    stop("Error: The input series must be an xts object.")
  }

  # Handle zero-break case explicitly
  if (is.null(breakpoints_dates) || length(breakpoints_dates) == 0) {
    cat("[create_dummy_variables] No breakpoints supplied. Returning empty list.\n")
    return(list())
  }

  time_index <- index(series)

  # Convert breakpoint dates to row indices in the series
  breakpoints <- match(breakpoints_dates, time_index)

  if (any(is.na(breakpoints))) {
    stop("Error: Some breakpoint dates are not found in the series index.")
  }

  n <- length(breakpoints)
  dummy_vars <- vector("list", n)

  partition_pts <- c(0, breakpoints, length(series))

  for (i in seq_len(n)) {
    dummy <- rep(0, length(series))
    dummy[(partition_pts[i + 1] + 1):partition_pts[i + 2]] <- 1
    dummy_vars[[i]] <- xts(dummy, order.by = time_index)
  }

  cat(sprintf("[create_dummy_variables] Created %d dummy variable(s) for %d break(s).\n", n, n))
  return(dummy_vars)
}

convert_monthly_to_daily <- function(monthly_series, trading_dates = NULL) {
  if (!inherits(monthly_series, "xts")) {
    stop("Error: The input series must be an xts object.")
  }

  monthly_dates <- index(monthly_series)
  values <- as.numeric(monthly_series)

  if (!is.null(trading_dates)) {
    trading_dates <- as.Date(trading_dates)
    trading_months <- floor_date(trading_dates, "month")
    monthly_df <- data.frame(
      month = floor_date(monthly_dates, "month"),
      value = values
    )
    matched_values <- monthly_df$value[match(trading_months, monthly_df$month)]

    # Handle NAs (trading dates outside monthly series range)
    if (any(is.na(matched_values))) {
      warning("Some trading dates have no matching monthly value. Setting to 0.")
      matched_values[is.na(matched_values)] <- 0
    }

    daily_series <- xts(matched_values, order.by = trading_dates)
    return(daily_series)
  }

  daily_data <- list()
  for (i in seq_along(monthly_dates)) {
    start_date <- as.Date(monthly_dates[i])
    end_date <- ceiling_date(start_date, "month") - 1
    daily_dates <- seq(start_date, end_date, by = "day")
    daily_values <- rep(values[i], length(daily_dates))
    daily_data[[i]] <- xts(daily_values, order.by = daily_dates)
  }

  daily_series <- do.call(rbind, daily_data)
  return(daily_series)
}

add_xts_matrix <- function(D_series, row, col, start_date, end_date) {
  XXX <- D_series[paste(start_date, end_date, sep = "/")]
  XXX_vec <- as.numeric(XXX)

  actual_len <- length(XXX_vec)

  if (actual_len != col) {
    stop(sprintf(
      "add_xts_matrix dimension mismatch: expected length %d (= number of trading days/cols), got %d.
         Dummy should have ONE value per trading day, which will be replicated across %d lag rows.
         Check that D_series has the same length as ncol(MV[[1]]).",
      col, actual_len, row
    ))
  }

  output_matrix <- matrix(rep(XXX_vec, each = row), nrow = row, ncol = col, byrow = FALSE)

  return(output_matrix)
}

# 2. MEAN MODEL ESTIMATION FUNCTION ==========================================================================================================================================
mean_fit <- function(returns, max_ar_order = 5, lb_lags = 10, alpha = 0.05,
                     force_ar_order = NULL, verbose = TRUE) {
  if (!"xts" %in% class(returns)) {
    stop("Input data must be an xts object.")
  }
  return_values <- as.numeric(coredata(returns))
  date_index <- index(returns)
  n <- length(return_values)
  valid_idx <- !is.na(return_values)
  rv_clean <- return_values[valid_idx]

  lb_raw <- Box.test(rv_clean, lag = lb_lags, type = "Ljung-Box")

  if (verbose) {
    cat("\n=== Mean Model Diagnostics ===\n")
    cat("Series:", if (!is.null(colnames(returns))) colnames(returns) else "unnamed", "\n")
    cat("Observations:", n, "\n")
    cat("\n[Step 1] Ljung-Box test on raw returns (lag =", lb_lags, "):\n")
    cat("  Statistic:", round(lb_raw$statistic, 4), "\n")
    cat("  p-value:  ", round(lb_raw$p.value, 4), "\n")
  }

  if (!is.null(force_ar_order) && force_ar_order > 0) {
    if (verbose) cat("  -> FORCED AR(", force_ar_order, ") by user.\n")
    best_order <- force_ar_order
    ar_fit_try <- tryCatch(
      arima(rv_clean, order = c(best_order, 0, 0), method = "ML"),
      error = function(e) NULL)
    if (is.null(ar_fit_try)) stop("Forced AR(", force_ar_order, ") failed to converge.")
    best_fit <- ar_fit_try
    model_type <- "AR"
    ar_order <- best_order
    ar_coefs <- coef(best_fit)[1:ar_order]
    mu_hat <- coef(best_fit)["intercept"]
    fitted_vals <- rv_clean - residuals(best_fit)
    residuals_full <- rep(NA, n)
    residuals_full[valid_idx] <- residuals(best_fit)
    residuals_est <- residuals_full

    std_error <- sqrt(best_fit$var.coef["intercept", "intercept"])
    t_stat <- mu_hat / std_error
    p_value <- 2 * (1 - pnorm(abs(t_stat)))

  } else if (lb_raw$p.value > alpha) {
    if (verbose) {
      cat("  -> No significant autocorrelation. Using constant mean.\n")
    }

    model_type <- "constant"
    ar_order <- 0
    mu_hat <- mean(rv_clean)
    residuals_est <- return_values - mu_hat 

    # Stats for constant mean
    std_error <- sd(rv_clean) / sqrt(length(rv_clean))
    t_stat <- mu_hat / std_error
    p_value <- 2 * (1 - pt(abs(t_stat), df = length(rv_clean) - 1))
    ar_coefs <- NULL

  } else {
    if (verbose) {
      cat("  -> Significant autocorrelation detected. Fitting AR(p) with BIC selection.\n")
    }
    bic_vals <- rep(NA, max_ar_order + 1)
    bic_vals[1] <- BIC(lm(rv_clean ~ 1)) 

    ar_fits <- list()
    for (p in 1:max_ar_order) {
      ar_fit_try <- tryCatch(
        arima(rv_clean, order = c(p, 0, 0), method = "ML"),
        error = function(e) NULL
      )
      if (!is.null(ar_fit_try)) {
        ar_fits[[p]] <- ar_fit_try
        bic_vals[p + 1] <- BIC(ar_fit_try)
      }
    }

    best_order <- which.min(bic_vals) - 1 

    if (verbose) {
      cat("  BIC values: AR(0):", round(bic_vals[1], 2))
      for (p in 1:max_ar_order) {
        if (!is.na(bic_vals[p + 1])) {
          cat(" AR(", p, "):", round(bic_vals[p + 1], 2), sep = "")
        }
      }
      cat("\n  -> Selected AR(", best_order, ")\n", sep = "")
    }

    if (best_order == 0) {
      model_type <- "constant"
      ar_order <- 0
      mu_hat <- mean(rv_clean)
      residuals_est <- return_values - mu_hat
      std_error <- sd(rv_clean) / sqrt(length(rv_clean))
      t_stat <- mu_hat / std_error
      p_value <- 2 * (1 - pt(abs(t_stat), df = length(rv_clean) - 1))
      ar_coefs <- NULL
    } else {
      model_type <- paste0("AR(", best_order, ")")
      ar_order <- best_order
      ar_fit <- ar_fits[[best_order]]

      mu_hat <- as.numeric(ar_fit$coef["intercept"])
      ar_coefs <- ar_fit$coef[1:best_order]

      resid_clean <- as.numeric(residuals(ar_fit))

      residuals_est <- rep(NA, n)
      residuals_est[valid_idx] <- resid_clean

      # Stats for intercept
      std_error <- sqrt(diag(vcov(ar_fit)))["intercept"]
      t_stat <- mu_hat / std_error
      p_value <- 2 * (1 - pnorm(abs(t_stat)))
    }
  }

  lb_resid <- Box.test(na.omit(residuals_est), lag = lb_lags, type = "Ljung-Box")

  if (verbose) {
    cat("\n[Step 3] Ljung-Box test on residuals (lag =", lb_lags, "):\n")
    cat("  Statistic:", round(lb_resid$statistic, 4), "\n")
    cat("  p-value:  ", round(lb_resid$p.value, 4), "\n")
    if (lb_resid$p.value > alpha) {
      cat("  -> Residuals are clean (no serial correlation).\n")
    } else {
      cat("  -> WARNING: Residuals still show serial correlation. Consider higher AR order.\n")
    }

    cat("\n[Summary]\n")
    cat("  Model type:", model_type, "\n")
    cat("  Intercept (mu_hat):", round(mu_hat, 6), "\n")
    cat("  SE:", round(std_error, 6), "\n")
    cat("  t-stat:", round(t_stat, 4), "\n")
    cat("  p-value:", round(p_value, 4), "\n")
    cat("==================================\n\n")
  }

  # Return results 
  residuals_xts <- xts(residuals_est, order.by = date_index)
  colnames(residuals_xts) <- "Residuals"

  return(list(
    estimated_mean = mu_hat,
    residuals = residuals_xts,
    standard_error = std_error,
    t_statistic = t_stat,
    p_value = p_value,
    model_type = model_type,
    ar_order = ar_order,
    ar_coefficients = ar_coefs,
    lb_test_raw = lb_raw,
    lb_test_residuals = lb_resid
  ))
}

# 3. GARCH-MIDAS ESTIMATION FUNCTIONS =============================================================================================================================================
GM_LL_gen <- function(param, daily_ret, mv_m, mv_m_m = NULL, K,
                       distribution = "norm", lag_fun = "Beta",
                       dummy = NULL, has_X = TRUE) {

  n_breaks <- if (is.null(dummy)) 0L else length(dummy)
  mu       <- param[1]
  alpha    <- param[2]
  beta     <- param[3]
  m        <- param[4]
  theta_RV <- param[5]
  w2_RV    <- param[6]
  w1       <- ifelse(lag_fun == "Beta", 1, 0)

  if (has_X) {
    theta_X  <- param[7]
    w2_X     <- param[8]
    base_idx <- 8L
  } else {
    base_idx <- 6L
  }

  if (n_breaks > 0) {
    theta_RV_b <- param[(base_idx + 1):(base_idx + n_breaks)]
    if (has_X) {
      theta_X_b <- param[(base_idx + n_breaks + 1):(base_idx + 2*n_breaks)]
    }
  }

  TT   <- length(daily_ret)
  g_it <- rep(1, TT)
  weight_fun <- ifelse(lag_fun == "Beta", beta_function, exp_almon)

  # MIDAS components
  betas_RV <- c(rev(weight_fun(1:(K+1), (K+1), w1, w2_RV))[2:(K+1)], 0)
  RV_comp  <- suppressWarnings(roll_sum(mv_m, c(K+1), weights = betas_RV))
  rv_sb <- 0
  if (n_breaks > 0) {
    for (j in 1:n_breaks) rv_sb <- rv_sb + theta_RV_b[j] * dummy[[j]]
  }

  if (has_X) {
    betas_X <- c(rev(weight_fun(1:(K+1), (K+1), w1, w2_X))[2:(K+1)], 0)
    XX_comp <- suppressWarnings(roll_sum(mv_m_m, c(K+1), weights = betas_X))

    x_sb <- 0
    if (n_breaks > 0) {
      for (j in 1:n_breaks) x_sb <- x_sb + theta_X_b[j] * dummy[[j]]
    }

    tau_d <- exp(pmin(pmax(m + (theta_RV + rv_sb) * RV_comp + (theta_X + x_sb) * XX_comp, -20), 20))
  } else {
    tau_d <- exp(pmin(pmax(m + (theta_RV + rv_sb) * RV_comp, -20), 20))
  }
  tau_d <- as.numeric(tau_d[K+1, ])

  step_1 <- (1 - alpha - beta) + alpha * (daily_ret - mu)^2 / tau_d
  for (i in 2:TT) {
    g_it[i] <- sum(step_1[i-1], beta * g_it[i-1], na.rm = TRUE)
  }

  ll <- as.numeric(stats::dnorm(daily_ret, mu, sqrt(g_it * tau_d), log = TRUE))
  return(ll)
}

GM_cond_vol_gen <- function(param, daily_ret, mv_m, mv_m_m = NULL, K,
                             lag_fun = "Beta", dummy = NULL, has_X = TRUE) {

  n_breaks <- if (is.null(dummy)) 0L else length(dummy)

  mu       <- param[1]
  alpha    <- param[2]; beta <- param[3]; m <- param[4]
  theta_RV <- param[5]; w2_RV <- param[6]
  w1       <- ifelse(lag_fun == "Beta", 1, 0)

  if (has_X) { theta_X <- param[7]; w2_X <- param[8]; base_idx <- 8L }
  else       { base_idx <- 6L }

  if (n_breaks > 0) {
    theta_RV_b <- param[(base_idx + 1):(base_idx + n_breaks)]
    if (has_X) theta_X_b <- param[(base_idx + n_breaks + 1):(base_idx + 2*n_breaks)]
  }

  TT   <- length(daily_ret); g_it <- rep(1, TT)
  weight_fun <- ifelse(lag_fun == "Beta", beta_function, exp_almon)
  betas_RV   <- c(rev(weight_fun(1:(K+1), (K+1), w1, w2_RV))[2:(K+1)], 0)
  RV_comp    <- suppressWarnings(roll_sum(mv_m, c(K+1), weights = betas_RV))

  rv_sb <- 0
  if (n_breaks > 0) for (j in 1:n_breaks) rv_sb <- rv_sb + theta_RV_b[j] * dummy[[j]]

  if (has_X) {
    betas_X <- c(rev(weight_fun(1:(K+1), (K+1), w1, w2_X))[2:(K+1)], 0)
    XX_comp <- suppressWarnings(roll_sum(mv_m_m, c(K+1), weights = betas_X))
    x_sb <- 0
    if (n_breaks > 0) for (j in 1:n_breaks) x_sb <- x_sb + theta_X_b[j] * dummy[[j]]
    tau_d <- exp(pmin(pmax(m + (theta_RV + rv_sb) * RV_comp + (theta_X + x_sb) * XX_comp, -20), 20))
  } else {
    tau_d <- exp(pmin(pmax(m + (theta_RV + rv_sb) * RV_comp, -20), 20))
  }
  tau_d <- as.numeric(tau_d[K+1, ])

  step_1 <- (1 - alpha - beta) + alpha * (daily_ret - mu)^2 / tau_d
  for (i in 2:TT) g_it[i] <- sum(step_1[i-1], beta * g_it[i-1], na.rm = TRUE)

  return(as.xts(sqrt(g_it * tau_d), index(daily_ret)))
}

GM_long_run_vol_gen <- function(param, daily_ret, mv_m, mv_m_m = NULL, K,
                                 lag_fun = "Beta", dummy = NULL, has_X = TRUE) {

  n_breaks <- if (is.null(dummy)) 0L else length(dummy)

  mu       <- param[1]
  alpha    <- param[2]; beta <- param[3]; m <- param[4]
  theta_RV <- param[5]; w2_RV <- param[6]
  w1       <- ifelse(lag_fun == "Beta", 1, 0)

  if (has_X) { theta_X <- param[7]; w2_X <- param[8]; base_idx <- 8L }
  else       { base_idx <- 6L }

  if (n_breaks > 0) {
    theta_RV_b <- param[(base_idx + 1):(base_idx + n_breaks)]
    if (has_X) theta_X_b <- param[(base_idx + n_breaks + 1):(base_idx + 2*n_breaks)]
  }

  TT <- length(daily_ret)
  weight_fun <- ifelse(lag_fun == "Beta", beta_function, exp_almon)
  betas_RV   <- c(rev(weight_fun(1:(K+1), (K+1), w1, w2_RV))[2:(K+1)], 0)
  RV_comp    <- suppressWarnings(roll_sum(mv_m, c(K+1), weights = betas_RV))

  rv_sb <- 0
  if (n_breaks > 0) for (j in 1:n_breaks) rv_sb <- rv_sb + theta_RV_b[j] * dummy[[j]]

  if (has_X) {
    betas_X <- c(rev(weight_fun(1:(K+1), (K+1), w1, w2_X))[2:(K+1)], 0)
    XX_comp <- suppressWarnings(roll_sum(mv_m_m, c(K+1), weights = betas_X))
    x_sb <- 0
    if (n_breaks > 0) for (j in 1:n_breaks) x_sb <- x_sb + theta_X_b[j] * dummy[[j]]
    tau_d <- exp(pmin(pmax(m + (theta_RV + rv_sb) * RV_comp + (theta_X + x_sb) * XX_comp, -20), 20))
  } else {
    tau_d <- exp(pmin(pmax(m + (theta_RV + rv_sb) * RV_comp, -20), 20))
  }
  tau_d <- as.numeric(tau_d[K+1, ])

  return(as.xts(sqrt(tau_d), index(daily_ret)))
}

ugmfit_gen <- function(daily_ret, mv_m, mv_m_2 = NULL, K,
                        lag_fun = "Beta", dummy = NULL, has_X = TRUE,
                        R = 800, seed = 1234) {

  n_breaks <- if (is.null(dummy)) 0L else length(dummy)
  n_params <- 6L + 2L*as.integer(has_X) + n_breaks*(1L + as.integer(has_X))

  # --- Build parameter names dynamically ---
  pnames <- c("mu", "alpha", "beta", "m", "theta_RV", "w2_RV")
  if (has_X) pnames <- c(pnames, "theta_X", "w2_X")
  if (n_breaks > 0) {
    pnames <- c(pnames, paste0("theta_RV_break_", 1:n_breaks))
    if (has_X) pnames <- c(pnames, paste0("theta_X_break_", 1:n_breaks))
  }

  cat(sprintf("\n[ugmfit_gen] has_X=%s, n_breaks=%d, n_params=%d\n",
              has_X, n_breaks, n_params))

  # --- Initial-value grid ---
  set.seed(seed)
  begin_val <- matrix(NA, nrow = R, ncol = n_params)
  colnames(begin_val) <- pnames

  begin_val[, 1] <- runif(R, -0.3, 0.5)    
  begin_val[, 2] <- runif(R, 0.001, 0.09)   
  begin_val[, 3] <- runif(R, 0.5, 0.85)     
  begin_val[, 4] <- runif(R, -5, 5)         
  begin_val[, 5] <- runif(R, -0.1, 0.1)    
  begin_val[, 6] <- runif(R, 1.05, 8.0)    

  if (has_X) {
    begin_val[, 7] <- runif(R, -5, 5)     
    begin_val[, 8] <- runif(R, 1.05, 8.0)   
  }

  base_idx <- if (has_X) 8L else 6L
  if (n_breaks > 0) {
    for (j in 1:n_breaks) {
      begin_val[, base_idx + j] <- runif(R, -0.1, 0.1)  
    }
    if (has_X) {
      for (j in 1:n_breaks) {
        begin_val[, base_idx + n_breaks + j] <- runif(R, -5, 5)   
      }
    }
  }

  r_t_est <- zoo::coredata(daily_ret)
  ll_args <- list(daily_ret = r_t_est, mv_m = mv_m, K = K,
                  distribution = "norm", lag_fun = lag_fun,
                  dummy = dummy, has_X = has_X)
  if (has_X) ll_args$mv_m_m <- mv_m_2

  which_row <- sapply(1:R, function(i) {
    sum(do.call(GM_LL_gen, c(list(param = begin_val[i, ]), ll_args)))
  })
  N_start <- max(3L, min(10L, 2L + n_breaks * 2L + as.integer(has_X)))
  top_idx <- order(which_row, decreasing = TRUE)[1:min(N_start, R)]

  ui_list <- list()
  ci_vec  <- c()

  row1 <- rep(0, n_params); row1[2] <- 1
  ui_list[[length(ui_list) + 1]] <- row1; ci_vec <- c(ci_vec, -1e-3)
  row2 <- rep(0, n_params); row2[3] <- 1
  ui_list[[length(ui_list) + 1]] <- row2; ci_vec <- c(ci_vec, -1e-3)
  row3 <- rep(0, n_params); row3[2] <- -1; row3[3] <- -1
  ui_list[[length(ui_list) + 1]] <- row3; ci_vec <- c(ci_vec, 0.999)
  row4 <- rep(0, n_params); row4[6] <- 1
  ui_list[[length(ui_list) + 1]] <- row4; ci_vec <- c(ci_vec, -1.001)
  
  if (has_X) {
    row5 <- rep(0, n_params); row5[8] <- 1
    ui_list[[length(ui_list) + 1]] <- row5; ci_vec <- c(ci_vec, -1.001)
  }

  row_ub1 <- rep(0, n_params); row_ub1[6] <- -1
  ui_list[[length(ui_list) + 1]] <- row_ub1; ci_vec <- c(ci_vec, 30)
  if (has_X) {
    row_ub2 <- rep(0, n_params); row_ub2[8] <- -1
    ui_list[[length(ui_list) + 1]] <- row_ub2; ci_vec <- c(ci_vec, 30)
  }
  ui <- do.call(rbind, ui_list)

  neg_ll_pen <- function(p) {
    ll_val <- sum(do.call(GM_LL_gen, c(list(param = p), ll_args)))
    if (!is.finite(ll_val)) return(1e10)
    pen <- 0
    if (p[2] < 1e-4) pen <- pen + 1e6 * (1e-4 - p[2])^2  
    if (p[3] < 1e-4) pen <- pen + 1e6 * (1e-4 - p[3])^2  
    if (p[2] + p[3] > 0.999) pen <- pen + 1e6 * (p[2]+p[3]-0.999)^2
    if (p[6] < 1.001) pen <- pen + 1e6 * (1.001 - p[6])^2
    if (p[6] > 30)    pen <- pen + 1e6 * (p[6] - 30)^2
    if (has_X && length(p) >= 8) {
      if (p[8] < 1.001) pen <- pen + 1e6 * (1.001 - p[8])^2
      if (p[8] > 30)    pen <- pen + 1e6 * (p[8] - 30)^2
    }
    -(ll_val - pen)
  }

  best_est <- NULL; best_ll <- -Inf
  for (s in seq_along(top_idx)) {
    sv <- begin_val[top_idx[s], ]
    nm_res <- tryCatch(optim(sv, fn = neg_ll_pen, method = "Nelder-Mead",
                             control = list(maxit = 3000)), error = function(e) NULL)
    if (!is.null(nm_res)) sv <- nm_res$par
    est_try <- tryCatch(suppressWarnings(maxLik::maxLik(
      logLik = GM_LL_gen, start = sv,
      daily_ret = r_t_est, mv_m = mv_m, mv_m_m = mv_m_2, K = K,
      distribution = "norm", lag_fun = lag_fun, dummy = dummy, has_X = has_X,
      constraints = list(ineqA = ui, ineqB = ci_vec),
      iterlim = 3000, method = "BFGS")), error = function(e) NULL)
    if (!is.null(est_try)) {
      ll_try <- as.numeric(logLik(est_try))
      if (ll_try > best_ll) { best_ll <- ll_try; best_est <- est_try }
      cat(sprintf("  [start %d/%d] LogL=%.4f%s\n", s, length(top_idx), ll_try,
                  if (ll_try >= best_ll) " <-- best" else ""))
    }
  }
  if (is.null(best_est)) stop("ugmfit_gen: all optimization attempts failed")
  est <- best_est
  est_coef <- coef(est)
  names(est_coef) <- pnames

  se_vals <- qmle_inference_engine(est, ll_func = GM_LL_gen, ll_args = ll_args,
                          daily_ret = r_t_est)
                  
  t_check <- abs(est_coef / se_vals)
  suspect <- which(!is.na(t_check) & is.finite(t_check) & t_check > 200)
  if (length(suspect) > 0) {
    cat("  [SE WARNING] Implausible t-stats (>200) for:", paste(pnames[suspect], collapse=", "),
        "\n        Marking these SEs as NA (Hessian likely unreliable at boundary).\n")
    se_vals[suspect] <- NA
  }

  mat_coef <- data.frame(
    Estimate     = round(est_coef, 6),
    `Std. Error` = round(se_vals, 6),
    `t value`    = round(est_coef / se_vals, 6),
    `Pr(>|t|)`   = round(2 * (1 - pnorm(abs(est_coef / se_vals))), 6),
    check.names = FALSE)
  rownames(mat_coef) <- pnames

  # --- Fitted volatility series ---
  vol_full <- GM_cond_vol_gen(est_coef, daily_ret, mv_m, mv_m_2, K,
                               lag_fun, dummy, has_X)
  lr_vol   <- GM_long_run_vol_gen(est_coef, daily_ret, mv_m, mv_m_2, K,
                                   lag_fun, dummy, has_X)

  loglik_val <- as.numeric(logLik(est))
  T_obs <- length(r_t_est)

  cat(sprintf("[ugmfit_gen] DONE: LogL=%.4f, AIC=%.4f, BIC=%.4f\n\n",
              loglik_val,
              -2*loglik_val + 2*n_params,
              -2*loglik_val + n_params*log(T_obs)))

  list(rob_coef_mat = mat_coef,
       est_vol_in_s = vol_full,
       est_lr_in_s  = lr_vol,
       loglik       = loglik_val,
       inf_criteria = c(AIC = -2*loglik_val + 2*n_params,
                        BIC = -2*loglik_val + n_params*log(T_obs)),
       n_params     = n_params,
       n_breaks     = n_breaks,
       has_X        = has_X,
       est_obj      = est)
}

# 4. DCC-MIDAS ESTIMATION FUNCTIONS ============================================================================================================================================
dccmidas_ll <- function (param, res, lag_fun = "Beta", N_c, K_c, dcc_mv = NULL)

{
  a <- param[1]
  b <- param[2]
  w1 <- ifelse(lag_fun == "Beta", 1, 0)
  w2 <- param[3]
  m <- param[4]
  theta <- param[5]
  theta_m <- param[6]
  w1_m <- ifelse(lag_fun == "Beta", 1, 0)
  w2_m <- param[7]
  Num_col <- dim(res)[1]
  TT <- dim(res)[3] 


  X_t <- array(0, dim = c(Num_col, Num_col, TT))
  for (tt in 1:TT) {
    X_t[, , tt] <- matrix(dcc_mv[tt], nrow = Num_col, ncol = Num_col)
  }

  C_t <- array(0, dim = c(Num_col, Num_col, TT))     
  V_t <- array(0, dim = c(Num_col, Num_col, TT))    
  Prod_eps_t <- array(0, dim = c(Num_col, Num_col, TT)) 

  S_init <- stats::cov(t(apply(res, 3L, c)))
  Q_t <- array(S_init, dim = c(Num_col, Num_col, TT))
  R_t <- array(diag(rep(1, Num_col)), dim = c(Num_col, Num_col, TT))

  log_det_R_t <- rep(0, TT)
  R_t_solved <- array(1, dim = c(Num_col, Num_col, TT))
  Eps_t_cross_prod <- rep(0, TT)
  Eps_t_R_t_Eps_t <- rep(0, TT)

  #Calculate Realized Correlation
  for (tt in (N_c + 1):TT) {
    V_t[, , tt] <- rowSums(res[, 1, tt:(tt - N_c)] * (res[, 1, tt:(tt - N_c)])) * diag(Num_col) 
    Prod_eps_t[, , tt] <- res[, , tt:(tt - N_c)] %*% t(res[, , tt:(tt - N_c)])
    V_t_0.5 <- Inv(sqrt(V_t[, , tt]))
    C_t[, , tt] <- V_t_0.5 %*% Prod_eps_t[, , tt] %*% V_t_0.5
  }

  weight_fun <- ifelse(lag_fun == "Beta", rumidas::beta_function, rumidas::exp_almon)
  betas <- c(rev(weight_fun(1:(K_c + 1), (K_c + 1), w1, w2))[2:(K_c +  1)], 0)
  weight_fun_m <- ifelse(lag_fun == "Beta", rumidas::beta_function, rumidas::exp_almon)
  betas_m <- c(rev(weight_fun_m(1:(K_c + 1), (K_c + 1), w1_m, w2_m))[2:(K_c +  1)], 0) #GEPU

  R_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  Z_t_bar <- array(1, dim = c(Num_col, Num_col, TT)) 
  matrix_id <- matrix(1:Num_col^2, ncol = Num_col)
  matrix_id_2 <- which(matrix_id == 1:Num_col^2, arr.ind = TRUE)

  C_t[, , 1:N_c] <- NA

  for (i in 1:nrow(matrix_id_2)) {
    C_component <- suppressWarnings(roll::roll_sum(C_t[matrix_id_2[i, 1], matrix_id_2[i, 2], ], c(K_c + 1), weights = betas))
    X_component <- suppressWarnings(roll::roll_sum(X_t[matrix_id_2[i, 1], matrix_id_2[i, 2], ], c(K_c + 1), weights = betas_m))
    Z_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ] <- m + theta*C_component + theta_m*(X_component)
    R_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ] <- (exp(2 * Z_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ]) - 1) /
      (exp(2 * Z_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ]) + 1)
  }

  burn_in <- N_c + K_c
  Z_t_bar[, , 1:burn_in] <- diag(Num_col)
  R_t_bar[, , 1:burn_in] <- diag(Num_col)

  for (tt in 1:dim(Z_t_bar)[3]) {
    diag(Z_t_bar[, , tt]) <- 1
    diag(R_t_bar[, , tt]) <- 1
  }

  ll <- rep(0, TT)


  for (tt in (burn_in + 1):TT) {
    Q_t[, , tt] <- (1 - a - b) * R_t_bar[, , tt] + a * res[, , tt - 1] %*% t(res[, , tt - 1]) + b * Q_t[, , tt - 1]
    Q_t_star <- Inv(sqrt(diag(diag(Q_t[, , tt]))))
    R_t[, , tt] <- Q_t_star %*% Q_t[, , tt] %*% Q_t_star

    log_det_R_t[tt] <- log(Det(R_t[, , tt]))
    R_t_solved[, , tt] <- Inv(R_t[, , tt])
    Eps_t_R_t_Eps_t[tt] <- rbind(res[, , tt]) %*% R_t_solved[, , tt] %*% cbind(res[, , tt])
  }

  ll <- -(log_det_R_t + Eps_t_R_t_Eps_t)
  return(ll)
}

# --- Log-likelihood with structural breaks ---
dccmidas_ll_sb <- function (param, res, lag_fun = "Beta", N_c, K_c, dcc_mv = NULL, dummy_dcc = NULL)
{
  a <- param[1]
  b <- param[2]
  w1 <- ifelse(lag_fun == "Beta", 1, 0)
  w2 <- param[3]
  m <- param[4]
  theta <- param[5]       
  theta_m <- param[6]     
  w1_m <- ifelse(lag_fun == "Beta", 1, 0)
  w2_m <- param[7]
  n_breaks <- if (is.null(dummy_dcc)) 0L else length(dummy_dcc)

  if (n_breaks > 0) {
    theta_b   <- param[(7 + 1):(7 + n_breaks)]
    theta_m_b <- param[(7 + n_breaks + 1):(7 + 2*n_breaks)]
  }

  Num_col <- dim(res)[1]; TT <- dim(res)[3]

  theta_sb   <- 0
  theta_m_sb <- 0
  if (n_breaks > 0) {
    for (j in 1:n_breaks) {
      theta_sb   <- theta_sb   + theta_b[j]   * dummy_dcc[[j]]
      theta_m_sb <- theta_m_sb + theta_m_b[j] * dummy_dcc[[j]]
    }
  }

  X_t <- array(0, dim = c(Num_col, Num_col, TT))
  for (tt in 1:TT) {
    X_t[, , tt] <- matrix(dcc_mv[tt], nrow = Num_col, ncol = Num_col)
  }
  C_t <- array(0, dim = c(Num_col, Num_col, TT))
  V_t <- array(0, dim = c(Num_col, Num_col, TT))
  Prod_eps_t <- array(0, dim = c(Num_col, Num_col, TT))

  S_init <- stats::cov(t(apply(res, 3L, c)))
  Q_t <- array(S_init, dim = c(Num_col, Num_col, TT))
  R_t <- array(diag(rep(1, Num_col)), dim = c(Num_col, Num_col, TT))

  log_det_R_t <- rep(0, TT)
  R_t_solved <- array(1, dim = c(Num_col, Num_col, TT))
  Eps_t_R_t_Eps_t <- rep(0, TT)

  for (tt in (N_c + 1):TT) {
    V_t[, , tt] <- rowSums(res[, 1, tt:(tt - N_c)] * (res[, 1, tt:(tt - N_c)])) * diag(Num_col)
    Prod_eps_t[, , tt] <- res[, , tt:(tt - N_c)] %*% t(res[, , tt:(tt - N_c)])
    V_t_0.5 <- Inv(sqrt(V_t[, , tt]))
    C_t[, , tt] <- V_t_0.5 %*% Prod_eps_t[, , tt] %*% V_t_0.5
  }

  weight_fun <- ifelse(lag_fun == "Beta", rumidas::beta_function, rumidas::exp_almon)
  betas <- c(rev(weight_fun(1:(K_c + 1), (K_c + 1), w1, w2))[2:(K_c + 1)], 0)
  weight_fun_m <- ifelse(lag_fun == "Beta", rumidas::beta_function, rumidas::exp_almon)
  betas_m <- c(rev(weight_fun_m(1:(K_c + 1), (K_c + 1), w1_m, w2_m))[2:(K_c + 1)], 0)

  R_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  Z_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  matrix_id <- matrix(1:Num_col^2, ncol = Num_col)
  matrix_id_2 <- which(matrix_id == 1:Num_col^2, arr.ind = TRUE)

  C_t[, , 1:N_c] <- NA

  for (i in 1:nrow(matrix_id_2)) {
    C_component <- suppressWarnings(roll::roll_sum(C_t[matrix_id_2[i, 1], matrix_id_2[i, 2], ], c(K_c + 1), weights = betas))
    X_component <- suppressWarnings(roll::roll_sum(X_t[matrix_id_2[i, 1], matrix_id_2[i, 2], ], c(K_c + 1), weights = betas_m))

    Z_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ] <- m +
      (theta + theta_sb) * C_component +
      (theta_m + theta_m_sb) * X_component

    R_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ] <- (exp(2 * Z_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ]) - 1) /
      (exp(2 * Z_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ]) + 1)
  }

  burn_in <- N_c + K_c
  Z_t_bar[, , 1:burn_in] <- diag(Num_col)
  R_t_bar[, , 1:burn_in] <- diag(Num_col)

  for (tt in 1:dim(Z_t_bar)[3]) {
    diag(Z_t_bar[, , tt]) <- 1
    diag(R_t_bar[, , tt]) <- 1
  }

  ll <- rep(0, TT)

  for (tt in (burn_in + 1):TT) {
    Q_t[, , tt] <- (1 - a - b) * R_t_bar[, , tt] + a * res[, , tt - 1] %*% t(res[, , tt - 1]) + b * Q_t[, , tt - 1]
    Q_t_star <- Inv(sqrt(diag(diag(Q_t[, , tt]))))
    R_t[, , tt] <- Q_t_star %*% Q_t[, , tt] %*% Q_t_star

    log_det_R_t[tt] <- log(Det(R_t[, , tt]))
    R_t_solved[, , tt] <- Inv(R_t[, , tt])
    Eps_t_R_t_Eps_t[tt] <- rbind(res[, , tt]) %*% R_t_solved[, , tt] %*% cbind(res[, , tt])
  }

  ll <- -(log_det_R_t + Eps_t_R_t_Eps_t)
  return(ll)
}


# --- Matrix extraction with structural breaks ---
dccmidas_mat_est_sb <- function (param, res, Dt, lag_fun = "Beta", N_c, K_c, dcc_mv = NULL, dummy_dcc = NULL)
{
  a <- param[1]; b <- param[2]
  w1 <- ifelse(lag_fun == "Beta", 1, 0); w2 <- param[3]
  m <- param[4]; theta <- param[5]; theta_m <- param[6]
  w1_m <- ifelse(lag_fun == "Beta", 1, 0); w2_m <- param[7]
  n_breaks <- if (is.null(dummy_dcc)) 0L else length(dummy_dcc)

  if (n_breaks > 0) {
    theta_b   <- param[(7 + 1):(7 + n_breaks)]
    theta_m_b <- param[(7 + n_breaks + 1):(7 + 2*n_breaks)]
  }

  Num_col <- dim(res)[1]; TT <- dim(res)[3]

  theta_sb   <- 0
  theta_m_sb <- 0
  if (n_breaks > 0) {
    for (j in 1:n_breaks) {
      theta_sb   <- theta_sb   + theta_b[j]   * dummy_dcc[[j]]
      theta_m_sb <- theta_m_sb + theta_m_b[j] * dummy_dcc[[j]]
    }
  }

  X_t <- array(0, dim = c(Num_col, Num_col, TT))
  for (tt in 1:TT) {
    X_t[, , tt] <- matrix(dcc_mv[tt], nrow = Num_col, ncol = Num_col)
  }

  C_t <- array(0, dim = c(Num_col, Num_col, TT))
  V_t <- array(0, dim = c(Num_col, Num_col, TT))
  Prod_eps_t <- array(0, dim = c(Num_col, Num_col, TT))
  S_init <- stats::cov(t(apply(res, 3L, c)))
  Q_t <- array(S_init, dim = c(Num_col, Num_col, TT))
  R_t <- array(diag(rep(1, Num_col)), dim = c(Num_col, Num_col, TT))
  S <- stats::cov(t(apply(res, 3L, c)))
  H_t <- array(S, dim = c(Num_col, Num_col, TT))

  for (tt in (N_c + 1):TT) {
    V_t[, , tt] <- rowSums(res[, 1, tt:(tt - N_c)] * (res[, 1, tt:(tt - N_c)])) * diag(Num_col)
    Prod_eps_t[, , tt] <- res[, , tt:(tt - N_c)] %*% t(res[, , tt:(tt - N_c)])
    V_t_0.5 <- Inv(sqrt(V_t[, , tt]))
    C_t[, , tt] <- V_t_0.5 %*% Prod_eps_t[, , tt] %*% V_t_0.5
  }

  weight_fun <- ifelse(lag_fun == "Beta", beta_function, exp_almon)
  betas <- c(rev(weight_fun(1:(K_c + 1), (K_c + 1), w1, w2))[2:(K_c + 1)], 0)
  weight_fun_m <- ifelse(lag_fun == "Beta", rumidas::beta_function, rumidas::exp_almon)
  betas_m <- c(rev(weight_fun_m(1:(K_c + 1), (K_c + 1), w1_m, w2_m))[2:(K_c + 1)], 0)

  Z_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  R_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  matrix_id <- matrix(1:Num_col^2, ncol = Num_col)
  matrix_id_2 <- which(matrix_id == 1:Num_col^2, arr.ind = TRUE)

  C_t[, , 1:N_c] <- NA

  for (i in 1:nrow(matrix_id_2)) {
    C_component <- suppressWarnings(roll::roll_sum(C_t[matrix_id_2[i, 1], matrix_id_2[i, 2], ], c(K_c + 1), weights = betas))
    X_component <- suppressWarnings(roll::roll_sum(X_t[matrix_id_2[i, 1], matrix_id_2[i, 2], ], c(K_c + 1), weights = betas_m))

    Z_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ] <- m +
      (theta + theta_sb) * C_component +
      (theta_m + theta_m_sb) * X_component

    R_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ] <- (exp(2 * Z_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ]) - 1) /
      (exp(2 * Z_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ]) + 1)
  }

  burn_in <- N_c + K_c
  Z_t_bar[, , 1:burn_in] <- diag(Num_col)
  R_t_bar[, , 1:burn_in] <- diag(Num_col)

  for (tt in 1:dim(Z_t_bar)[3]) {
    diag(Z_t_bar[, , tt]) <- 1
    diag(R_t_bar[, , tt]) <- 1
  }

  for (tt in (burn_in + 1):TT) {
    Q_t[, , tt] <- (1 - a - b) * R_t_bar[, , tt] + a * res[, , tt - 1] %*% t(res[, , tt - 1]) + b * Q_t[, , tt - 1]
    Q_t_star <- Inv(sqrt(diag(diag(Q_t[, , tt]))))
    R_t[, , tt] <- Q_t_star %*% Q_t[, , tt] %*% Q_t_star
    H_t[, , tt] <- Dt[, , tt] %*% R_t[, , tt] %*% Dt[, , tt]
  }

  results <- list(H_t = H_t, R_t = R_t, R_t_bar = R_t_bar, C_t = C_t, Z_t = Z_t_bar,
                  X_t = X_t, Q_t = Q_t, D_t = Dt, res = res, C_component = C_component, X_component = X_component)
  return(results)
}

dccmidas_mat_est <- function (param, res, Dt, lag_fun = "Beta", N_c, K_c, dcc_mv = NULL)
{
  a <- param[1]
  b <- param[2]
  w1 <- ifelse(lag_fun == "Beta", 1, 0)
  w2 <- param[3]
  m <- param[4]
  theta <- param[5]
  theta_m <- param[6]
  w1_m <- ifelse(lag_fun == "Beta", 1, 0)
  w2_m <- param[7]
  Num_col <- dim(res)[1] 
  TT <- dim(res)[3] 

  X_t <- array(0, dim = c(Num_col, Num_col, TT))
  for (tt in 1:TT) {
    X_t[, , tt] <- matrix(dcc_mv[tt], nrow = Num_col, ncol = Num_col)
  } 


  C_t <- array(0, dim = c(Num_col, Num_col, TT))  
  V_t <- array(0, dim = c(Num_col, Num_col, TT))  

  Prod_eps_t <- array(0, dim = c(Num_col, Num_col, TT))  
  S_init <- stats::cov(t(apply(res, 3L, c)))
  Q_t <- array(S_init, dim = c(Num_col, Num_col, TT))
  R_t <- array(diag(rep(1, Num_col)), dim = c(Num_col, Num_col, TT)) 

  log_det_R_t <- rep(0, TT)  
  R_t_solved <- array(1, dim = c(Num_col, Num_col, TT))  
  Eps_t_cross_prod <- rep(0, TT)  
  Eps_t_R_t_Eps_t <- rep(0, TT)  
  S <- stats::cov(t(apply(res, 3L, c)))  
  H_t <- array(S, dim = c(Num_col, Num_col, TT)) 

  # Compute rolling covariance and normalize
  for (tt in (N_c + 1):TT) {
    V_t[, , tt] <- rowSums(res[, 1, tt:(tt - N_c)] * (res[, 1, tt:(tt - N_c)])) * diag(Num_col)
    Prod_eps_t[, , tt] <- res[, , tt:(tt - N_c)] %*% t(res[, , tt:(tt - N_c)])
    V_t_0.5 <- Inv(sqrt(V_t[, , tt]))
    C_t[, , tt] <- V_t_0.5 %*% Prod_eps_t[, , tt] %*% V_t_0.5
  }

  weight_fun <- ifelse(lag_fun == "Beta", beta_function, exp_almon)
  betas <- c(rev(weight_fun(1:(K_c + 1), (K_c + 1), w1, w2))[2:(K_c + 1)], 0)
  weight_fun_m <- ifelse(lag_fun == "Beta", rumidas::beta_function, rumidas::exp_almon)
  betas_m <- c(rev(weight_fun_m(1:(K_c + 1), (K_c + 1), w1_m, w2_m))[2:(K_c +  1)], 0)

  Z_t_bar <- array(1, dim = c(Num_col, Num_col, TT)) 
  R_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  matrix_id <- matrix(1:Num_col^2, ncol = Num_col)
  matrix_id_2 <- which(matrix_id == 1:Num_col^2, arr.ind = TRUE)

  C_t[, , 1:N_c] <- NA

  for (i in 1:nrow(matrix_id_2)) {
    C_component <- suppressWarnings(roll::roll_sum(C_t[matrix_id_2[i, 1], matrix_id_2[i, 2], ], c(K_c + 1), weights = betas))
    X_component <- suppressWarnings(roll::roll_sum(X_t[matrix_id_2[i, 1], matrix_id_2[i, 2], ], c(K_c + 1), weights = betas_m))
    Z_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ] <- m + theta*C_component + theta_m*(X_component)
    R_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ] <- (exp(2 * Z_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ]) - 1) /
      (exp(2 * Z_t_bar[matrix_id_2[i, 1], matrix_id_2[i, 2], ]) + 1)
  }

  burn_in <- N_c + K_c
  Z_t_bar[, , 1:burn_in] <- diag(Num_col)
  R_t_bar[, , 1:burn_in] <- diag(Num_col)

  for (tt in 1:dim(Z_t_bar)[3]) {
    diag(Z_t_bar[, , tt]) <- 1
    diag(R_t_bar[, , tt]) <- 1
  }

  for (tt in (burn_in + 1):TT) {
    Q_t[, , tt] <- (1 - a - b) * R_t_bar[, , tt] + a * res[, , tt - 1] %*% t(res[, , tt - 1]) + b * Q_t[, , tt - 1]
    Q_t_star <- Inv(sqrt(diag(diag(Q_t[, , tt]))))
    R_t[, , tt] <- Q_t_star %*% Q_t[, , tt] %*% Q_t_star 
    H_t[, , tt] <- Dt[, , tt] %*% R_t[, , tt] %*% Dt[, , tt] 
  }

  results <- list(H_t = H_t, R_t = R_t, R_t_bar = R_t_bar, C_t = C_t, Z_t = Z_t_bar, X_t = X_t, Q_t = Q_t, D_t = Dt, res = res, C_component = C_component, X_component = X_component)
  return(results)
}


#  RC-ONLY DCC-MIDAS FUNCTIONS (NO GEPU)                                             

# ---------- Log-likelihood: RC only, WITHOUT SB ----------

dccmidas_ll_rc_no_sb <- function(param, res, lag_fun = "Beta", N_c, K_c) {
  a <- param[1]; b <- param[2]
  w1 <- ifelse(lag_fun == "Beta", 1, 0); w2 <- param[3]
  m <- param[4]; theta <- param[5]

  Num_col <- dim(res)[1]; TT <- dim(res)[3]

  C_t <- array(0, dim = c(Num_col, Num_col, TT))
  V_t <- array(0, dim = c(Num_col, Num_col, TT))
  Prod_eps_t <- array(0, dim = c(Num_col, Num_col, TT))
  S_init <- stats::cov(t(apply(res, 3L, c)))
  Q_t <- array(S_init, dim = c(Num_col, Num_col, TT))
  R_t <- array(diag(rep(1, Num_col)), dim = c(Num_col, Num_col, TT))
  log_det_R_t <- rep(0, TT)
  R_t_solved <- array(1, dim = c(Num_col, Num_col, TT))
  Eps_t_R_t_Eps_t <- rep(0, TT)

  for (tt in (N_c+1):TT) {
    V_t[,,tt] <- rowSums(res[,1,tt:(tt-N_c)] * res[,1,tt:(tt-N_c)]) * diag(Num_col)
    Prod_eps_t[,,tt] <- res[,,tt:(tt-N_c)] %*% t(res[,,tt:(tt-N_c)])
    V_t_0.5 <- Inv(sqrt(V_t[,,tt]))
    C_t[,,tt] <- V_t_0.5 %*% Prod_eps_t[,,tt] %*% V_t_0.5
  }

  weight_fun <- ifelse(lag_fun == "Beta", rumidas::beta_function, rumidas::exp_almon)
  betas <- c(rev(weight_fun(1:(K_c+1), (K_c+1), w1, w2))[2:(K_c+1)], 0)

  R_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  Z_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  matrix_id <- matrix(1:Num_col^2, ncol = Num_col)
  matrix_id_2 <- which(matrix_id == 1:Num_col^2, arr.ind = TRUE)

  C_t[, , 1:N_c] <- NA

  for (i in 1:nrow(matrix_id_2)) {
    C_component <- suppressWarnings(roll::roll_sum(
      C_t[matrix_id_2[i,1], matrix_id_2[i,2], ], c(K_c+1), weights = betas))
    Z_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ] <- m + theta * C_component
    R_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ] <-
      (exp(2 * Z_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ]) - 1) /
      (exp(2 * Z_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ]) + 1)
  }

  burn_in <- N_c + K_c
  Z_t_bar[, , 1:burn_in] <- diag(Num_col)
  R_t_bar[, , 1:burn_in] <- diag(Num_col)
  for (tt in 1:dim(Z_t_bar)[3]) {
    diag(Z_t_bar[,,tt]) <- 1; diag(R_t_bar[,,tt]) <- 1
  }

  for (tt in (burn_in + 1):TT) {
    Q_t[,,tt] <- (1-a-b)*R_t_bar[,,tt] + a*res[,,tt-1]%*%t(res[,,tt-1]) + b*Q_t[,,tt-1]
    Q_t_star <- Inv(sqrt(diag(diag(Q_t[,,tt]))))
    R_t[,,tt] <- Q_t_star %*% Q_t[,,tt] %*% Q_t_star
    log_det_R_t[tt] <- log(Det(R_t[,,tt]))
    R_t_solved[,,tt] <- Inv(R_t[,,tt])
    Eps_t_R_t_Eps_t[tt] <- rbind(res[,,tt]) %*% R_t_solved[,,tt] %*% cbind(res[,,tt])
  }

  return(-(log_det_R_t + Eps_t_R_t_Eps_t))
}

# ---------- Matrix extraction: RC only, WITHOUT SB ----------
dccmidas_mat_est_rc_no_sb <- function(param, res, Dt, lag_fun = "Beta", N_c, K_c) {
  a <- param[1]; b <- param[2]
  w1 <- ifelse(lag_fun == "Beta", 1, 0); w2 <- param[3]
  m <- param[4]; theta <- param[5]

  Num_col <- dim(res)[1]; TT <- dim(res)[3]
  C_t <- array(0, dim = c(Num_col, Num_col, TT))
  V_t <- array(0, dim = c(Num_col, Num_col, TT))
  Prod_eps_t <- array(0, dim = c(Num_col, Num_col, TT))
  S_init <- stats::cov(t(apply(res, 3L, c)))
  Q_t <- array(S_init, dim = c(Num_col, Num_col, TT))
  R_t <- array(diag(rep(1, Num_col)), dim = c(Num_col, Num_col, TT))
  S <- stats::cov(t(apply(res, 3L, c)))
  H_t <- array(S, dim = c(Num_col, Num_col, TT))

  for (tt in (N_c+1):TT) {
    V_t[,,tt] <- rowSums(res[,1,tt:(tt-N_c)] * res[,1,tt:(tt-N_c)]) * diag(Num_col)
    Prod_eps_t[,,tt] <- res[,,tt:(tt-N_c)] %*% t(res[,,tt:(tt-N_c)])
    V_t_0.5 <- Inv(sqrt(V_t[,,tt]))
    C_t[,,tt] <- V_t_0.5 %*% Prod_eps_t[,,tt] %*% V_t_0.5
  }

  weight_fun <- ifelse(lag_fun == "Beta", rumidas::beta_function, rumidas::exp_almon)
  betas <- c(rev(weight_fun(1:(K_c+1), (K_c+1), w1, w2))[2:(K_c+1)], 0)

  Z_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  R_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  matrix_id <- matrix(1:Num_col^2, ncol = Num_col)
  matrix_id_2 <- which(matrix_id == 1:Num_col^2, arr.ind = TRUE)

  C_t[, , 1:N_c] <- NA

  for (i in 1:nrow(matrix_id_2)) {
    C_component <- suppressWarnings(roll::roll_sum(
      C_t[matrix_id_2[i,1], matrix_id_2[i,2], ], c(K_c+1), weights = betas))
    Z_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ] <- m + theta * C_component
    R_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ] <-
      (exp(2 * Z_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ]) - 1) /
      (exp(2 * Z_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ]) + 1)
  }

  burn_in <- N_c + K_c
  Z_t_bar[, , 1:burn_in] <- diag(Num_col)
  R_t_bar[, , 1:burn_in] <- diag(Num_col)
  for (tt in 1:dim(Z_t_bar)[3]) {
    diag(Z_t_bar[,,tt]) <- 1; diag(R_t_bar[,,tt]) <- 1
  }

  for (tt in (burn_in + 1):TT) {
    Q_t[,,tt] <- (1-a-b)*R_t_bar[,,tt] + a*res[,,tt-1]%*%t(res[,,tt-1]) + b*Q_t[,,tt-1]
    Q_t_star <- Inv(sqrt(diag(diag(Q_t[,,tt]))))
    R_t[,,tt] <- Q_t_star %*% Q_t[,,tt] %*% Q_t_star
    H_t[,,tt] <- Dt[,,tt] %*% R_t[,,tt] %*% Dt[,,tt]
  }

  return(list(H_t=H_t, R_t=R_t, R_t_bar=R_t_bar, C_t=C_t, Z_t=Z_t_bar, Q_t=Q_t, D_t=Dt, res=res))
}

# ---------- Log-likelihood: RC only, WITH SB ----------

dccmidas_ll_rc_sb <- function(param, res, lag_fun = "Beta", N_c, K_c, dummy_dcc = NULL) {

  n_breaks <- if (is.null(dummy_dcc)) 0L else length(dummy_dcc)
  a <- param[1]; b <- param[2]
  w1 <- ifelse(lag_fun == "Beta", 1, 0); w2 <- param[3]
  m <- param[4]; theta <- param[5]

  if (n_breaks > 0) {
    theta_b <- param[(5 + 1):(5 + n_breaks)]   
  }

  Num_col <- dim(res)[1]; TT <- dim(res)[3]

  theta_sb <- 0
  if (n_breaks > 0) {
    for (j in 1:n_breaks) theta_sb <- theta_sb + theta_b[j] * dummy_dcc[[j]]
  }

  C_t <- array(0, dim = c(Num_col, Num_col, TT))
  V_t <- array(0, dim = c(Num_col, Num_col, TT))
  Prod_eps_t <- array(0, dim = c(Num_col, Num_col, TT))
  S_init <- stats::cov(t(apply(res, 3L, c)))
  Q_t <- array(S_init, dim = c(Num_col, Num_col, TT))
  R_t <- array(diag(rep(1, Num_col)), dim = c(Num_col, Num_col, TT))
  log_det_R_t <- rep(0, TT)
  R_t_solved <- array(1, dim = c(Num_col, Num_col, TT))
  Eps_t_R_t_Eps_t <- rep(0, TT)

  for (tt in (N_c+1):TT) {
    V_t[,,tt] <- rowSums(res[,1,tt:(tt-N_c)] * res[,1,tt:(tt-N_c)]) * diag(Num_col)
    Prod_eps_t[,,tt] <- res[,,tt:(tt-N_c)] %*% t(res[,,tt:(tt-N_c)])
    V_t_0.5 <- Inv(sqrt(V_t[,,tt]))
    C_t[,,tt] <- V_t_0.5 %*% Prod_eps_t[,,tt] %*% V_t_0.5
  }

  weight_fun <- ifelse(lag_fun == "Beta", rumidas::beta_function, rumidas::exp_almon)
  betas <- c(rev(weight_fun(1:(K_c+1), (K_c+1), w1, w2))[2:(K_c+1)], 0)

  R_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  Z_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  matrix_id <- matrix(1:Num_col^2, ncol = Num_col)
  matrix_id_2 <- which(matrix_id == 1:Num_col^2, arr.ind = TRUE)

  C_t[, , 1:N_c] <- NA

  for (i in 1:nrow(matrix_id_2)) {
    C_component <- suppressWarnings(roll::roll_sum(
      C_t[matrix_id_2[i,1], matrix_id_2[i,2], ], c(K_c+1), weights = betas))
    Z_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ] <- m + (theta + theta_sb) * C_component
    R_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ] <-
      (exp(2 * Z_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ]) - 1) /
      (exp(2 * Z_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ]) + 1)
  }

  burn_in <- N_c + K_c
  Z_t_bar[,,1:burn_in] <- diag(Num_col)
  R_t_bar[,,1:burn_in] <- diag(Num_col)
  for (tt in 1:dim(Z_t_bar)[3]) {
    diag(Z_t_bar[,,tt]) <- 1; diag(R_t_bar[,,tt]) <- 1
  }

  for (tt in (burn_in + 1):TT) {
    Q_t[,,tt] <- (1-a-b)*R_t_bar[,,tt] + a*res[,,tt-1]%*%t(res[,,tt-1]) + b*Q_t[,,tt-1]
    Q_t_star <- Inv(sqrt(diag(diag(Q_t[,,tt]))))
    R_t[,,tt] <- Q_t_star %*% Q_t[,,tt] %*% Q_t_star
    log_det_R_t[tt] <- log(Det(R_t[,,tt]))
    R_t_solved[,,tt] <- Inv(R_t[,,tt])
    Eps_t_R_t_Eps_t[tt] <- rbind(res[,,tt]) %*% R_t_solved[,,tt] %*% cbind(res[,,tt])
  }

  return(-(log_det_R_t + Eps_t_R_t_Eps_t))
}

# ---------- Matrix extraction: RC only, WITH SB ----------
dccmidas_mat_est_rc_sb <- function(param, res, Dt, lag_fun = "Beta", N_c, K_c, dummy_dcc = NULL) {

  n_breaks <- if (is.null(dummy_dcc)) 0L else length(dummy_dcc)

  a <- param[1]; b <- param[2]
  w1 <- ifelse(lag_fun == "Beta", 1, 0); w2 <- param[3]
  m <- param[4]; theta <- param[5]

  if (n_breaks > 0) {
    theta_b <- param[(5 + 1):(5 + n_breaks)]
  }

  Num_col <- dim(res)[1]; TT <- dim(res)[3]

  theta_sb <- 0
  if (n_breaks > 0) {
    for (j in 1:n_breaks) theta_sb <- theta_sb + theta_b[j] * dummy_dcc[[j]]
  }

  C_t <- array(0, dim = c(Num_col, Num_col, TT))
  V_t <- array(0, dim = c(Num_col, Num_col, TT))
  Prod_eps_t <- array(0, dim = c(Num_col, Num_col, TT))
  S_init <- stats::cov(t(apply(res, 3L, c)))
  Q_t <- array(S_init, dim = c(Num_col, Num_col, TT))
  R_t <- array(diag(rep(1, Num_col)), dim = c(Num_col, Num_col, TT))
  S <- stats::cov(t(apply(res, 3L, c)))
  H_t <- array(S, dim = c(Num_col, Num_col, TT))

  for (tt in (N_c+1):TT) {
    V_t[,,tt] <- rowSums(res[,1,tt:(tt-N_c)] * res[,1,tt:(tt-N_c)]) * diag(Num_col)
    Prod_eps_t[,,tt] <- res[,,tt:(tt-N_c)] %*% t(res[,,tt:(tt-N_c)])
    V_t_0.5 <- Inv(sqrt(V_t[,,tt]))
    C_t[,,tt] <- V_t_0.5 %*% Prod_eps_t[,,tt] %*% V_t_0.5
  }

  weight_fun <- ifelse(lag_fun == "Beta", rumidas::beta_function, rumidas::exp_almon)
  betas <- c(rev(weight_fun(1:(K_c+1), (K_c+1), w1, w2))[2:(K_c+1)], 0)

  Z_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  R_t_bar <- array(1, dim = c(Num_col, Num_col, TT))
  matrix_id <- matrix(1:Num_col^2, ncol = Num_col)
  matrix_id_2 <- which(matrix_id == 1:Num_col^2, arr.ind = TRUE)

  C_t[, , 1:N_c] <- NA

  for (i in 1:nrow(matrix_id_2)) {
    C_component <- suppressWarnings(roll::roll_sum(
      C_t[matrix_id_2[i,1], matrix_id_2[i,2], ], c(K_c+1), weights = betas))
    Z_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ] <- m + (theta + theta_sb) * C_component
    R_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ] <-
      (exp(2 * Z_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ]) - 1) /
      (exp(2 * Z_t_bar[matrix_id_2[i,1], matrix_id_2[i,2], ]) + 1)
  }

  burn_in <- N_c + K_c
  Z_t_bar[,,1:burn_in] <- diag(Num_col)
  R_t_bar[,,1:burn_in] <- diag(Num_col)
  for (tt in 1:dim(Z_t_bar)[3]) {
    diag(Z_t_bar[,,tt]) <- 1; diag(R_t_bar[,,tt]) <- 1
  }

  for (tt in (burn_in + 1):TT) {
    Q_t[,,tt] <- (1-a-b)*R_t_bar[,,tt] + a*res[,,tt-1]%*%t(res[,,tt-1]) + b*Q_t[,,tt-1]
    Q_t_star <- Inv(sqrt(diag(diag(Q_t[,,tt]))))
    R_t[,,tt] <- Q_t_star %*% Q_t[,,tt] %*% Q_t_star
    H_t[,,tt] <- Dt[,,tt] %*% R_t[,,tt] %*% Dt[,,tt]
  }

  return(list(H_t=H_t, R_t=R_t, R_t_bar=R_t_bar, C_t=C_t, Z_t=Z_t_bar, Q_t=Q_t, D_t=Dt, res=res))
}

dcc_fit_modify <- function (r_t,
                            univ_model = "GM_noskew", distribution = "norm",
                            MV = NULL, MV_M = NULL, dcc_mv = NULL, K = NULL,
                            corr_model = "DCCMIDAS", lag_fun = "Beta",
                            N_c = NULL, K_c = NULL, out_of_sample = NULL,
                            dummy = NULL)

{
  cond_r_t <- class(r_t)   
  if (cond_r_t != "list") {
    stop(cat("#Warning:\n Parameter 'r_t' must be a list object. Please provide it in the correct form \n"))
  }
  if ((univ_model %in% c("GM_noskew", "GM_skew", "DAGM_noskew", "DAGM_skew")) & (!inherits(MV, "list") | length(MV) != length(r_t) | is.null(K))) {
    stop(cat("#Warning:\n If you want to estimate a GARCH-MIDAS model, please provide parameters MV and K in the correct form \n"))
  }
  if ((corr_model == "DCCMIDAS" | corr_model == "ADCCMIDAS") &
      (is.null(N_c) | is.null(K_c))) {
    stop(cat("#Warning:\n If you want to estimate a DCC-MIDAS model, please provide parameters N_c and K_c \n"))
  }

  # ============================================================================
  Num_assets <- length(r_t)
  len <- rep(NA, Num_assets)

  for (i in 1:Num_assets) {
    len[i] <- length(r_t[[i]])
  }
  if (stats::var(len) != 0) {
    db <- do.call(xts::merge.xts, r_t)
    db <- db[stats::complete.cases(db), ]
  }
  else {
    db <- do.call(xts::merge.xts, r_t)
  }
  for (i in 1:Num_assets) {
    colnames(db)[i] <- colnames(r_t[[i]])
  }
  TT <- nrow(db)
  if (missing(out_of_sample)) {
    sd_m <- matrix(NA, ncol = Num_assets, nrow = TT)
    sd_m_lr <- matrix(NA, ncol = Num_assets, nrow = TT)
    estim_period <- range(stats::time(db))
    days_period <- stats::time(db)
    est_obs <- nrow(db)
  }
  else {
    sd_m <- matrix(NA, ncol = Num_assets, nrow = TT - out_of_sample)
    db_oos <- db[(TT - out_of_sample + 1):TT, ]
    sd_m_oos <- matrix(NA, ncol = Num_assets, nrow = out_of_sample)
    estim_period <- range(stats::time(db[1:(TT - out_of_sample),
    ]))
    days_period <- stats::time(db[1:(TT - out_of_sample),
    ])
    est_obs <- nrow(db[1:(TT - out_of_sample), ])
  }

  est_details <- list()
  likelihood_at_max <- list()
  start <- Sys.time()


  # ============================================================================

  garch_logL <- list()

  if (univ_model == "GM_noskew") {
    if (missing(out_of_sample)) {
      for (i in 1:Num_assets) {
        u_est <- ugmfit_gen(model = "GM", skew = "NO",
                             distribution = "norm", db[, i], mv_m = MV[[i]], mv_m_2 = MV_M[[i]],
                             K = K, lag_fun = lag_fun,
                             dummy = dummy)

        est_details[[i]] <- u_est$rob_coef_mat   # The robust coefficient matrix (rob_coef_mat) for asset
        likelihood_at_max[[i]] <- u_est$loglik   # The log-likelihood value at the estimated parameters' maximum
        sd_m[, i] <- u_est$est_vol_in_s          
        sd_m_lr[, i] <- u_est$est_lr_in_s
        garch_logL[[i]] <- u_est$inf_criteria
      }
    }
    else {
      for (i in 1:Num_assets) {
        u_est <- rumidas::ugmfit(model = "GM", skew = "NO",
                                 distribution = distribution, db[, i], MV[[i]],
                                 K = K, lag_fun = lag_fun, out_of_sample = out_of_sample)
        est_details[[i]] <- u_est$rob_coef_mat
        likelihood_at_max[[i]] <- u_est$loglik
        sd_m[, i] <- u_est$est_vol_in_s
        sd_m_oos[, i] <- zoo::coredata(u_est$est_vol_oos)
      }
    }
  }

  #=============================================================================


  cat("First step: completed \n")
  if (missing(out_of_sample)) {
    D_t <- array(0, dim = c(Num_assets, Num_assets, TT))
    eps_t <- array(0, dim = c(Num_assets, 1, TT))
    db_a <- array(NA, dim = c(Num_assets, 1, TT))
    db_no_xts <- zoo::coredata(db)
    for (tt in 1:TT) {
      db_a[, , tt] <- t(db_no_xts[tt, ])
      diag(D_t[, , tt]) <- sd_m[tt, ]
      eps_t[, , tt] <- Inv(D_t[, , tt]) %*% db_a[, , tt]
    }
  }
  else {
    D_t <- array(0, dim = c(Num_assets, Num_assets, TT -
                              out_of_sample))
    eps_t <- array(0, dim = c(Num_assets, 1, TT - out_of_sample))
    db_a <- array(NA, dim = c(Num_assets, 1, TT - out_of_sample))
    db_no_xts <- zoo::coredata(db[1:(TT - out_of_sample),
    ])
    for (tt in 1:(TT - out_of_sample)) {
      db_a[, , tt] <- t(db_no_xts[tt, ])
      diag(D_t[, , tt]) <- sd_m[tt, ]
      eps_t[, , tt] <- Inv(D_t[, , tt]) %*% db_a[, , tt]
    }
    D_t_oos <- array(0, dim = c(Num_assets, Num_assets, out_of_sample))
    eps_t_oos <- array(0, dim = c(Num_assets, 1, out_of_sample))
    db_a_oos <- array(NA, dim = c(Num_assets, 1, out_of_sample))
    db_no_xts_oos <- zoo::coredata(db_oos)
    for (tt in 1:out_of_sample) {
      db_a_oos[, , tt] <- t(db_no_xts_oos[tt, ])
      diag(D_t_oos[, , tt]) <- sd_m_oos[tt, ]
      eps_t_oos[, , tt] <- Inv(D_t_oos[, , tt]) %*% db_a_oos[,
                                                             , tt]
    }
  }

  #---------------------------------------------------

  if (corr_model == "DCCMIDAS" & lag_fun == "Beta") {
    R = 1000
    start_val <- begin_val <- ui <- ci <- NULL
    begin_val <- matrix(NA, nrow = R, ncol = 7)
    colnames(begin_val) <- c("alpha", "beta", "w2", "m", "theta", "theta_m", "w2_m")
    begin_val[, 1] <- stats::runif(R, min = 0.001, max = 0.095)
    begin_val[, 2] <- stats::runif(R, min = 0.3, max = 0.8)
    begin_val[, 3] <- runif(R, 1.05, 8.0)
    begin_val[, 4] <- stats::runif(R, min = -5, max = 5)
    begin_val[, 5] <- stats::runif(R, min = -1, max = 1)
    begin_val[, 6] <- stats::runif(R, min = -1, max = 1)
    begin_val[, 7] <- runif(R, 1.05, 8.0)

    which_row <- rep(NA, R)
    for (i in 1:R) {
      which_row[i] <- sum(dccmidas_ll(begin_val[i, ], res = eps_t, lag_fun = "Beta", N_c = N_c, K_c = K_c, dcc_mv = dcc_mv))
    }
    start_val <- begin_val[which.max(which_row), ]

    ui <- rbind(c(1, 0, 0, 0, 0, 0, 0),
                c(0, 1, 0, 0, 0, 0, 0),
                c(-1, -1, 0, 0, 0, 0, 0),
                c(0, 0, 1, 0, 0, 0, 0),
                c(0, 0, 0, 0, 0, 0, 1)
    )
    ci <- c(-1e-04, -0.001, 0.999, -1.001, -1.001)

    m_est <- maxLik(
      logLik = dccmidas_ll,
      start = start_val, res = eps_t,
      lag_fun = lag_fun, N_c = N_c, K_c = K_c, dcc_mv = dcc_mv,
      constraints = list(ineqA = ui, ineqB = ci),
      iterlim = 5000,
      method = "BFGS"
    )
  }

  #=============================================================================

  est_coef <- stats::coef(m_est)
  N_coef <- length(est_coef)
  mat_coef <- data.frame(rep(NA, N_coef), rep(NA, N_coef),
                         rep(NA, N_coef), rep(NA, N_coef))
  colnames(mat_coef) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(mat_coef) <- names(est_coef)
  mat_coef[, 1] <- round(est_coef, 6)
  se_dcc <- qmle_inference_engine(m_est,
                         ll_func = dccmidas_ll,
                         ll_args = list(res = eps_t, lag_fun = lag_fun, N_c = N_c, K_c = K_c, dcc_mv = dcc_mv))
  mat_coef[, 2] <- round(se_dcc, 6)
  mat_coef[, 3] <- round(est_coef / se_dcc, 6)
  mat_coef[, 4] <- round(
    apply(rbind(est_coef / se_dcc), 1,
          function(x) 2 * (1 - stats::pnorm(abs(x)))), 6)

  if (corr_model == "DCCMIDAS") {
    dcc_mat_est_fin <- dccmidas_mat_est(est_coef, eps_t, D_t, lag_fun = lag_fun, N_c = N_c, K_c = K_c, dcc_mv = dcc_mv)
    if (!missing(out_of_sample)) {
      dcc_mat_est_fin_oos <- dccmidas_mat_est(est_coef, eps_t_oos, D_t_oos, lag_fun = lag_fun, N_c = N_c, K_c = K_c, dcc_mv = dcc_mv)
    }
    else {
      dcc_mat_est_fin_oos <- list(NA, NA, NA)
    }
  }

  #=============================================================================

  end <- Sys.time() - start
  if (corr_model == "DCCMIDAS") {
    fin_res <- list(assets = colnames(db), model = univ_model,
                    est_univ_model = est_details, corr_coef_mat = mat_coef, garch_midas = est_details,
                    mult_model = corr_model, obs = est_obs, period = estim_period,
                    H_t = dcc_mat_est_fin[[1]], R_t = dcc_mat_est_fin[[2]], R_t_bar = dcc_mat_est_fin[[3]],
                    C_t = dcc_mat_est_fin[[4]], Z_t = dcc_mat_est_fin[[5]], X_t = dcc_mat_est_fin[[6]], Q_t = dcc_mat_est_fin[[7]],
                    D_t = dcc_mat_est_fin[[8]], res = dcc_mat_est_fin[[9]], C_component = dcc_mat_est_fin[[10]], X_component = dcc_mat_est_fin[[11]],
                    garch_midas_est_vol = sd_m, garch_midas_est_lvol = sd_m_lr, resid = db_a, garch_logL = garch_logL, likelihood_at_max = likelihood_at_max,
                    H_t_oos = dcc_mat_est_fin_oos[[1]], R_t_oos = dcc_mat_est_fin_oos[[2]], R_t_bar_oos = dcc_mat_est_fin_oos[[3]],
                    est_time = end, Days = days_period, llk = stats::logLik(m_est))
  }

  class(fin_res) <- c("dccmidas")
  return(fin_res)
  print.dccmidas(fin_res)
}

# 5. DATA FORMATTING FUNCTIONS ============================================================================================================================================

to_df_multiple <- function(...) {
  matrix_names <- as.character(match.call())[-1]
  matrices <- list(...)

  df_list <- mapply(function(mat, mat_name) {
    dims <- dim(mat)

    if (is.null(dims)) {
      return(data.frame(setNames(list(mat), mat_name)))
    }

    else if (length(dims) == 2) {
      temp_list <- setNames(as.data.frame(mat), paste0(mat_name, "_", seq_len(ncol(mat))))
      return(temp_list)
    }

    else if (length(dims) == 3) {
      temp_list <- list()
      for (i in seq_len(dims[1])) {
        for (j in seq_len(dims[2])) {
          col_name <- paste0(mat_name, "_", i, "_", j)
          temp_list[[col_name]] <- mat[i, j, ]
        }
      }
      return(as.data.frame(temp_list))
    }

    else {
      stop("Error: Only vectors, 2D, and 3D matrices are supported.")
    }
  }, matrices, matrix_names, SIMPLIFY = FALSE)

  return(do.call(cbind, df_list))
}
set_date_index <- function(df, start_date=start_date, end_date=end_date) {
  # Load necessary libraries
  library(dplyr)
  library(tibble)

  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)

  date_sequence <- seq(start_date, end_date, by = "day")

  if (nrow(df) != length(date_sequence)) {
    stop("Error: The number of rows in the data frame does not match the length of the date sequence.")
  }

  df$Date <- date_sequence

  df <- df %>% column_to_rownames(var = "Date")

  return(df)
}

# CONFIG 

INDEX_TICKER <- "msci_world"    # <<<<<<< CHANGE THIS for each pair global economic policy uncertainty

K_vol  <- 24   
K_corr <- 24    
N_c    <- 22   

start_date_is <- as.Date("2016-10-01")
end_date_is   <- as.Date("2025-11-30")

# 6. LOAD DATA FUNCTION ============================================================================================================================================

load_data <- function(start_date, end_date, index_ticker) {
  data <- read_excel("D:/Thesis/Data/BTC-INDICES/Aligned_Trading_Dates_All_Assets.xlsx")
  data <- data[order(as.Date(data$Date)), ]
  dates <- as.Date(data$Date)

  log_ret_btc <- as.numeric(data$Ln_rt_btc)*100
  log_ret_idx <- as.numeric(data[[paste0("Ln_rt_", index_ticker)]])*100

  btc_xts <- xts(log_ret_btc, order.by = dates)
  idx_xts <- xts(log_ret_idx, order.by = dates)
  colnames(btc_xts) <- "BTC"
  colnames(idx_xts) <- toupper(index_ticker)

  db_m <- merge.xts(btc_xts, idx_xts)
  db_m_filtered <- window(db_m, start = start_date, end = end_date)

  data_gepu <- read_excel("D:/Thesis/Data/BTC-INDICES/Data_Monthly_GEPU.xlsx", sheet = "GEPU")
  data_gepu <- data_gepu[order(as.Date(data_gepu$Date)), ]
  GEPU_xts <- xts(as.numeric(data_gepu$GEPU_ppp),
                  order.by = as.Date(data_gepu$Date))

  GEPU_diff <- diff(log(GEPU_xts))
  GEPU_diff <- GEPU_diff[!is.na(GEPU_diff)]
  colnames(GEPU_diff) <- "DLn_GEPU"
  cat("GEPU log-diff range:", round(min(GEPU_diff), 4), "to",
      round(max(GEPU_diff), 4), "\n")

  # --- Monthly RV ---
  data_rv <- read_excel("D:/Thesis/Data/BTC-INDICES/Data_Monthly_GEPU.xlsx", sheet = "RV")
  data_rv <- data_rv[order(as.Date(data_rv$Date)), ]
  dates_rv <- as.Date(data_rv$Date)

  RV_scale <- 100^2   # = 10000
  RV_btc_xts <- xts(as.numeric(data_rv$RV_Ln_rt_btc) * RV_scale, order.by = dates_rv)
  RV_idx_xts <- xts(as.numeric(data_rv[[paste0("RV_Ln_rt_", index_ticker)]]) * RV_scale,
                    order.by = dates_rv)
  cat(sprintf("RV scaled by %d to match percentage returns.\n", RV_scale))
  cat(sprintf("  BTC RV mean: %.4f  |  %s RV mean: %.4f\n",
              mean(RV_btc_xts, na.rm=TRUE), toupper(index_ticker), mean(RV_idx_xts, na.rm=TRUE)))

  # --- Transform to MIDAS matrix format ---
  r_t <- list(db_m_filtered[, 1], db_m_filtered[, 2])
  mv_m_btc <- mv_into_mat(r_t[[1]], RV_btc_xts, K = K_vol, "monthly")
  mv_m_idx <- mv_into_mat(r_t[[2]], RV_idx_xts, K = K_vol, "monthly")
  mv_m_m   <- mv_into_mat(r_t[[1]], GEPU_diff, K = K_vol, "monthly")

  MV   <- list(mv_m_btc, mv_m_idx)
  MV_M <- list(mv_m_m, mv_m_m)
  mv_x <- mv_into_mat(r_t[[1]], GEPU_diff, K = 0, "monthly")

  return(list(
    r_t = r_t, MV = MV, MV_M = MV_M, mv_x = mv_x,
    X = as.numeric(GEPU_diff), X_1 = GEPU_diff,
    X_level = as.numeric(GEPU_xts), X_level_xts = GEPU_xts
  ))
}

# 7. HEDGING EVALUATION FUNCTION ============================================================================================================================================

hedging_performance_analysis <- function(df) {
  if (!exists("df")) {
    stop("Error: The dataframe does not exist.")
  }
  
  df$hedge_ratio_2 <- df$H_t_1_2 / df$H_t_1_1
  df$hedge_ratio_1 <- df$H_t_1_2 / df$H_t_2_2

  df$ret_hedged_2 <- df$resid_2_1 - (df$hedge_ratio_2 * df$resid_1_1)
  df$ret_hedged_1 <- df$resid_1_1 - (df$hedge_ratio_1 * df$resid_2_1)

  compute_stats <- function(column_name) {
    var_value <- var(df[[column_name]], na.rm = TRUE) 
    hpr_value <- sum(df[[column_name]], na.rm = TRUE) 
    return(c(Variance = var_value, Return = hpr_value))
  }

  stats_resid_1_1 <- compute_stats("resid_1_1")
  stats_ret_hedged_1 <- compute_stats("ret_hedged_1")
  stats_resid_2_1 <- compute_stats("resid_2_1")
  stats_ret_hedged_2 <- compute_stats("ret_hedged_2")

  comparison_df <- data.frame(
    Metric = c("Variance", "Return"),
    Resid_1_1 = stats_resid_1_1,
    Ret_Hedged_1 = stats_ret_hedged_1,
    Resid_2_1 = stats_resid_2_1,
    Ret_Hedged_2 = stats_ret_hedged_2
  )

  # Extract variance values
  var_resid_1_1 <- as.numeric(comparison_df[comparison_df$Metric == "Variance", "Resid_1_1"])
  var_ret_hedged_1 <- as.numeric(comparison_df[comparison_df$Metric == "Variance", "Ret_Hedged_1"])
  var_resid_2_1 <- as.numeric(comparison_df[comparison_df$Metric == "Variance", "Resid_2_1"])
  var_ret_hedged_2 <- as.numeric(comparison_df[comparison_df$Metric == "Variance", "Ret_Hedged_2"])

  # Extract HPR values
  hpr_resid_1_1 <- as.numeric(comparison_df[comparison_df$Metric == "Return", "Resid_1_1"])
  hpr_ret_hedged_1 <- as.numeric(comparison_df[comparison_df$Metric == "Return", "Ret_Hedged_1"])
  hpr_resid_2_1 <- as.numeric(comparison_df[comparison_df$Metric == "Return", "Resid_2_1"])
  hpr_ret_hedged_2 <- as.numeric(comparison_df[comparison_df$Metric == "Return", "Ret_Hedged_2"])

  # Compute Hedging Effectiveness (HE)
  HE_1 <- (var_resid_1_1 - var_ret_hedged_1) / var_resid_1_1
  HE_2 <- (var_resid_2_1 - var_ret_hedged_2) / var_resid_2_1

  # Compute Sharpe Ratios
  sharpe_resid_1_1 <- hpr_resid_1_1 / sqrt(var_resid_1_1)
  sharpe_ret_hedged_1 <- hpr_ret_hedged_1 / sqrt(var_ret_hedged_1)
  sharpe_resid_2_1 <- hpr_resid_2_1 / sqrt(var_resid_2_1)
  sharpe_ret_hedged_2 <- hpr_ret_hedged_2 / sqrt(var_ret_hedged_2)

  # Create a dataframe for Sharpe Ratio comparison
  sharpe_comparison_df <- data.frame(
    Metric = c("Sharpe Ratio"),
    Resid_1_1 = sharpe_resid_1_1,
    Ret_Hedged_1 = sharpe_ret_hedged_1,
    Resid_2_1 = sharpe_resid_2_1,
    Ret_Hedged_2 = sharpe_ret_hedged_2
  )

  return(list(
    Comparison_Stats = comparison_df,
    Hedging_Effectiveness = c(HE_1 = HE_1, HE_2 = HE_2),
    Sharpe_Comparison = sharpe_comparison_df
  ))
}

build_eps_Dt <- function(sd_btc, sd_idx, db) {
  TT <- nrow(db)
  D_t   <- array(0, dim = c(2, 2, TT))
  eps_t <- array(0, dim = c(2, 1, TT))
  db_no_xts <- zoo::coredata(db)

  # Coerce to plain vectors
  sd_b <- as.numeric(sd_btc); sd_i <- as.numeric(sd_idx)

  for (tt in 1:TT) {
    diag(D_t[,,tt]) <- c(sd_b[tt], sd_i[tt])
    eps_t[,,tt] <- solve(D_t[,,tt]) %*% t(db_no_xts[tt, , drop = FALSE])
  }
  list(eps_t = eps_t, D_t = D_t, TT = TT)
}

# ============================================================================
estimate_dcc_gen <- function(ll_fun, mat_fun, eps_t, D_t, TT,
                              N_c, K_c, R_dcc = 1500,
                              has_X = FALSE, dummy_dcc = NULL,
                              dcc_mv = NULL, seed = 1234,
                              begin_val_override = NULL) {

  set.seed(seed)
  n_breaks <- if (is.null(dummy_dcc)) 0L else length(dummy_dcc)
  n_p <- 5L + 2L*as.integer(has_X) + n_breaks*(1L + as.integer(has_X))

  # --- Parameter names ---
  pnames <- c("a", "b", "w2_RC", "m", "theta_RC")
  if (has_X) pnames <- c(pnames, "theta_X", "w2_X")
  if (n_breaks > 0) {
    pnames <- c(pnames, paste0("theta_RC_break_", 1:n_breaks))
    if (has_X) pnames <- c(pnames, paste0("theta_X_break_", 1:n_breaks))
  }

  cat(sprintf("\n[estimate_dcc_gen] has_X=%s, n_breaks=%d, n_params=%d\n",
              has_X, n_breaks, n_p))

  # --- Initial-value grid ---
  begin_val <- matrix(NA, nrow = R_dcc, ncol = n_p)
  colnames(begin_val) <- pnames
  begin_val[, 1] <- runif(R_dcc, 0.001, 0.095) 
  begin_val[, 2] <- runif(R_dcc, 0.5, 0.98)      
  begin_val[, 3] <- runif(R_dcc, 1.05, 8.0)       
  begin_val[, 4] <- runif(R_dcc, -3, 3)          
  begin_val[, 5] <- runif(R_dcc, -5, 5)        
  if (has_X) {
    begin_val[, 6] <- runif(R_dcc, -10, 10)  
    begin_val[, 7] <- runif(R_dcc, 1.05, 8.0)  
  }
  base_idx <- if (has_X) 7L else 5L
  if (n_breaks > 0) {
    for (j in 1:n_breaks)
      begin_val[, base_idx + j] <- runif(R_dcc, -5, 5)
    if (has_X) {
      for (j in 1:n_breaks)
        begin_val[, base_idx + n_breaks + j] <- runif(R_dcc, -10, 10)
    }
  }

  # --- Build LL args ---
  base_args <- list(res = eps_t, lag_fun = "Beta", N_c = N_c, K_c = K_c)
  if (!is.null(dcc_mv))    base_args$dcc_mv    <- dcc_mv
  if (!is.null(dummy_dcc)) base_args$dummy_dcc <- dummy_dcc

  # --- Rank grid starting values ---
  which_row <- sapply(1:R_dcc, function(i) {
    args_i <- c(list(param = begin_val[i, ]), base_args)
    sum(do.call(ll_fun, args_i))
  })

  if (!is.null(begin_val_override)) {
    n_override <- length(begin_val_override)
    if (n_override <= n_p) {
      warm_row <- rep(0, n_p)
      warm_row[1:n_override] <- begin_val_override
      n_warm <- min(5L, R_dcc)
      worst_idx <- order(which_row, decreasing = FALSE)[1:n_warm]
      for (wi in seq_along(worst_idx)) {
        perturb <- warm_row * runif(n_p, 0.9, 1.1)
        perturb[1] <- max(0.001, min(0.095, perturb[1]))
        perturb[2] <- max(0.5, min(0.95, perturb[2]))
        perturb[3] <- max(1.05, min(29, perturb[3]))
        if (has_X && n_p >= 7) perturb[7] <- max(1.05, min(29, perturb[7]))
        begin_val[worst_idx[wi], ] <- perturb
      }
      begin_val[worst_idx[1], ] <- warm_row
      which_row <- sapply(1:R_dcc, function(i) {
        args_i <- c(list(param = begin_val[i, ]), base_args)
        sum(do.call(ll_fun, args_i))
      })
      cat("  [warm-start] Injected from simpler model.\n")
    }
  }

  N_start <- max(3L, min(10L, 2L + n_breaks * 2L + as.integer(has_X) * 2L))
  top_idx <- order(which_row, decreasing = TRUE)[1:min(N_start, R_dcc)]

  ui_list <- list()
  ci_vec  <- c()

  row1 <- rep(0, n_p); row1[1] <- 1
  ui_list[[length(ui_list)+1]] <- row1; ci_vec <- c(ci_vec, -1e-3)

  row2 <- rep(0, n_p); row2[2] <- 1
  ui_list[[length(ui_list)+1]] <- row2; ci_vec <- c(ci_vec, -0.5)
  row3 <- rep(0, n_p); row3[1] <- -1; row3[2] <- -1
  ui_list[[length(ui_list)+1]] <- row3; ci_vec <- c(ci_vec, 0.999)
  row4 <- rep(0, n_p); row4[3] <- 1
  ui_list[[length(ui_list)+1]] <- row4; ci_vec <- c(ci_vec, -1.001)
  if (has_X) {
    row5 <- rep(0, n_p); row5[7] <- 1
    ui_list[[length(ui_list)+1]] <- row5; ci_vec <- c(ci_vec, -1.001)
  }

  row_t1 <- rep(0, n_p); row_t1[5] <- 1
  ui_list[[length(ui_list)+1]] <- row_t1; ci_vec <- c(ci_vec, 20)  
  row_t2 <- rep(0, n_p); row_t2[5] <- -1
  ui_list[[length(ui_list)+1]] <- row_t2; ci_vec <- c(ci_vec, 20)   
  if (has_X) {

    row_t3 <- rep(0, n_p); row_t3[6] <- 1
    ui_list[[length(ui_list)+1]] <- row_t3; ci_vec <- c(ci_vec, 20)
    row_t4 <- rep(0, n_p); row_t4[6] <- -1
    ui_list[[length(ui_list)+1]] <- row_t4; ci_vec <- c(ci_vec, 20)
  }

  row_ub1 <- rep(0, n_p); row_ub1[3] <- -1
  ui_list[[length(ui_list)+1]] <- row_ub1; ci_vec <- c(ci_vec, 30)
  if (has_X) {
    row_ub2 <- rep(0, n_p); row_ub2[7] <- -1
    ui_list[[length(ui_list)+1]] <- row_ub2; ci_vec <- c(ci_vec, 30)
  }
  ui <- do.call(rbind, ui_list)

  neg_ll_pen_dcc <- function(p) {
    ll_val <- sum(do.call(ll_fun, c(list(param = p), base_args)))
    if (!is.finite(ll_val)) return(1e10)
    pen <- 0
    if (p[1] < 1e-4)  pen <- pen + 1e6 * (1e-4 - p[1])^2
    if (p[2] < 0.5)   pen <- pen + 1e6 * (0.5 - p[2])^2     
    if (p[1] + p[2] > 0.999) pen <- pen + 1e6 * (p[1]+p[2]-0.999)^2
    if (p[3] < 1.001) pen <- pen + 1e6 * (1.001 - p[3])^2
    if (p[3] > 30)    pen <- pen + 1e6 * (p[3] - 30)^2

    if (abs(p[5]) > 20) pen <- pen + 1e6 * (abs(p[5]) - 20)^2
    if (has_X && length(p) >= 7) {
      if (p[7] < 1.001) pen <- pen + 1e6 * (1.001 - p[7])^2
      if (p[7] > 30)    pen <- pen + 1e6 * (p[7] - 30)^2

      if (abs(p[6]) > 20) pen <- pen + 1e6 * (abs(p[6]) - 20)^2
    }
    -(ll_val - pen)
  }

  best_est <- NULL; best_ll <- -Inf
  for (s in seq_along(top_idx)) {
    sv <- begin_val[top_idx[s], ]
    nm_res <- tryCatch(optim(sv, fn = neg_ll_pen_dcc, method = "Nelder-Mead",
                             control = list(maxit = 5000)), error = function(e) NULL)
    if (!is.null(nm_res)) sv <- nm_res$par
    est_try <- tryCatch(suppressWarnings(do.call(maxLik, c(
      list(logLik = ll_fun, start = sv), base_args,
      list(constraints = list(ineqA = ui, ineqB = ci_vec),
           iterlim = 8000, method = "BFGS")))), error = function(e) NULL)
    if (!is.null(est_try)) {
      ll_try <- as.numeric(logLik(est_try))
      if (ll_try > best_ll) { best_ll <- ll_try; best_est <- est_try }
      cat(sprintf("  [start %d/%d] LogL=%.4f%s\n", s, length(top_idx), ll_try,
                  if (ll_try >= best_ll) " <-- best" else ""))
    }
  }
  if (is.null(best_est)) stop("estimate_dcc_gen: all optimization attempts failed")
  m_est <- best_est

  # --- Coefficient table ---
  est_coef <- coef(m_est)
  names(est_coef) <- pnames
  se <- qmle_inference_engine(m_est, ll_func = ll_fun, ll_args = base_args)

  t_check <- abs(est_coef / se)
  suspect <- which(!is.na(t_check) & is.finite(t_check) & t_check > 200)
  if (length(suspect) > 0) {
    cat("  [SE WARNING] Implausible t-stats (>200) for:", paste(pnames[suspect], collapse=", "),
        "\n        Marking these SEs as NA (Hessian likely unreliable at boundary).\n")
    se[suspect] <- NA
  }

  mat_coef <- data.frame(
    Estimate     = round(est_coef, 6),
    `Std. Error` = round(se, 6),
    `t value`    = round(est_coef / se, 6),
    `Pr(>|t|)`   = round(2 * (1 - pnorm(abs(est_coef / se))), 6),
    check.names = FALSE)
  rownames(mat_coef) <- pnames

  # --- Extract matrices ---
  mat_args <- list(param = est_coef, res = eps_t, Dt = D_t,
                   lag_fun = "Beta", N_c = N_c, K_c = K_c)
  if (!is.null(dcc_mv))    mat_args$dcc_mv    <- dcc_mv
  if (!is.null(dummy_dcc)) mat_args$dummy_dcc <- dummy_dcc
  dcc_matrices <- do.call(mat_fun, mat_args)

  llk <- as.numeric(logLik(m_est))
  list(coef_mat = mat_coef, matrices = dcc_matrices,
       llk = llk, AIC = -2*llk + 2*n_p, BIC = -2*llk + n_p*log(TT),
       n_params = n_p, n_breaks = n_breaks, has_X = has_X, m_est = m_est)
}

# B. =========================================================================== ------------ MAIN EXECUTION  ------------ =========================================================================  

# 1. LOAD DATA ===========================================================================================================================================
load_data_is <- load_data(start_date_is, end_date_is, INDEX_TICKER)
r_t       <- load_data_is$r_t
MV        <- load_data_is$MV
MV_M      <- load_data_is$MV_M
mv_x      <- load_data_is$mv_x
macro     <- load_data_is$X
macro_xts <- load_data_is$X_1

cat("\n========== PAIR: BTC -", toupper(INDEX_TICKER), "==========\n")

macro_level_xts <- load_data_is$X_level_xts

# 2. STRUCTURAL BREAK DETECTION ==========================================================================================================================
break_test <- bai_perron_test(
  series        = macro_level_xts,   
  max_breaks    = NULL,               
  selection     = "BIC",               # BIC-optimal selection
  min_supF_pval = 0.05,                # require supF significance
  start_date    = start_date_is,
  end_date      = end_date_is)

breakpoints_vec <- break_test$breakpoints_vector  
n_breaks        <- break_test$n_breaks             

dummy_list <- create_dummy_variables(
  breakpoints_dates = breakpoints_vec,
  series            = macro_level_xts)

cat(sprintf("\n>>> Detected %d structural break(s). Downstream SB models will use %d dummy(ies).\n\n",
            n_breaks, length(dummy_list)))

cat(">>> mu estimated jointly inside GARCH-MIDAS likelihood (no pre-demeaning)\n")

r_t_  <- list(r_t[[1]], r_t[[2]])

db_temp <- do.call(xts::merge.xts, r_t_)
db_temp <- db_temp[stats::complete.cases(db_temp), ]
trading_dates_final <- index(db_temp)

if (nrow(db_temp) != ncol(MV[[1]])) {
  cat("NOTE: Realigning MV due to NA rows in raw returns...\n")
  data_rv <- read_excel("D:/Thesis/Data/BTC-INDICES/Data_Monthly_GEPU.xlsx", sheet = "RV")
  data_rv <- data_rv[order(as.Date(data_rv$Date)), ]
  dates_rv <- as.Date(data_rv$Date)
  RV_scale <- 100^2
  RV_btc_xts <- xts(as.numeric(data_rv$RV_Ln_rt_btc) * RV_scale, order.by = dates_rv)
  RV_idx_xts <- xts(as.numeric(data_rv[[paste0("RV_Ln_rt_", INDEX_TICKER)]]) * RV_scale, order.by = dates_rv)

  data_gepu <- read_excel("D:/Thesis/Data/BTC-INDICES/Data_Monthly_GEPU.xlsx", sheet = "GEPU")
  data_gepu <- data_gepu[order(as.Date(data_gepu$Date)), ]
  GEPU_xts_level <- xts(as.numeric(data_gepu$GEPU_ppp), order.by = as.Date(data_gepu$Date))
  GEPU_diff_aligned <- diff(log(GEPU_xts_level))
  GEPU_diff_aligned <- GEPU_diff_aligned[!is.na(GEPU_diff_aligned)]

  MV   <- list(mv_into_mat(db_temp[,1], RV_btc_xts, K=K_vol, "monthly"),
               mv_into_mat(db_temp[,2], RV_idx_xts, K=K_vol, "monthly"))
  MV_M <- list(mv_into_mat(db_temp[,1], GEPU_diff_aligned, K=K_vol, "monthly"),
               mv_into_mat(db_temp[,1], GEPU_diff_aligned, K=K_vol, "monthly"))
  mv_x <- mv_into_mat(db_temp[,1], GEPU_diff_aligned, K=0, "monthly")
}

if (length(dummy_list) == 0) {
  dummy_input     <- NULL  
  dummy_dcc_input <- NULL   
  cat(">>> No structural breaks: SB tables will be identical to no-SB tables.\n\n")
} else {
  dummy_input <- list()
  for (i in seq_along(dummy_list)) {
    D_month <- dummy_list[[i]]
    D_daily <- convert_monthly_to_daily(D_month, trading_dates = trading_dates_final)
    D_final <- add_xts_matrix(D_daily, dim(MV[[1]])[1], dim(MV[[1]])[2],
                              start_date_is, end_date_is)
    dummy_input[[i]] <- D_final
  }

  dummy_dcc_input <- list()
  for (i in seq_along(dummy_list)) {
    D_month <- dummy_list[[i]]
    D_daily <- convert_monthly_to_daily(D_month, trading_dates = trading_dates_final)
    dummy_dcc_input[[i]] <- as.numeric(D_daily)
  }

  cat(sprintf(">>> Built %d daily dummy series for GARCH-MIDAS and DCC-MIDAS.\n\n",
              length(dummy_input)))
}

# 3. PREPARE DATA FOR GARCH-MIDAS AND DCC-MIDAS ===========================================================================================================================

# GARCH-MIDAS ESTIMATION
                  
cat("\n", strrep("=", 70), "\n")
cat(" GARCH-MIDAS ESTIMATION: 8 TABLES (gen wrapper)\n")
cat(strrep("=", 70), "\n")

cat("\n>>> [Tables 1-2] GARCH-MIDAS RV-only, no SB\n")
garch_rv_no_sb <- list()
for (i in 1:2) {
  garch_rv_no_sb[[i]] <- ugmfit_gen(
    daily_ret = r_t_[[i]], mv_m = MV[[i]],
    K = K_vol, has_X = FALSE, dummy = NULL, R = 800)
}

cat("\n>>> [Tables 3-4] GARCH-MIDAS RV-only, with SB\n")
garch_rv_sb <- list()
for (i in 1:2) {
  garch_rv_sb[[i]] <- ugmfit_gen(
    daily_ret = r_t_[[i]], mv_m = MV[[i]],
    K = K_vol, has_X = FALSE, dummy = dummy_input, R = 800)
}

cat("\n>>> [Tables 5-6] GARCH-MIDAS-X (RV+GEPU), no SB\n")
garch_rvx_no_sb <- list()
for (i in 1:2) {
  garch_rvx_no_sb[[i]] <- ugmfit_gen(
    daily_ret = r_t_[[i]], mv_m = MV[[i]], mv_m_2 = MV_M[[i]],
    K = K_vol, has_X = TRUE, dummy = NULL, R = 800)
}

cat("\n>>> [Tables 7-8] GARCH-MIDAS-X (RV+GEPU), with SB\n")
garch_rvx_sb <- list()
for (i in 1:2) {
  garch_rvx_sb[[i]] <- ugmfit_gen(
    daily_ret = r_t_[[i]], mv_m = MV[[i]], mv_m_2 = MV_M[[i]],
    K = K_vol, has_X = TRUE, dummy = dummy_input, R = 800)
}

# DCC-MIDAS ESTIMATION                                                     
# All DCC models use eps_t from the GARCH-MIDAS-X with SB (full model estimation)     

cat("\n" , strrep("=", 70), "\n")
cat(" DCC-MIDAS ESTIMATION: 4 TABLES\n")
cat(strrep("=", 70), "\n")

db <- do.call(xts::merge.xts, r_t_)
db <- db[stats::complete.cases(db), ]
TT <- nrow(db)

# --- 
ed_table9  <- build_eps_Dt(garch_rv_no_sb[[1]]$est_vol_in_s,
                           garch_rv_no_sb[[2]]$est_vol_in_s, db)
ed_table10 <- build_eps_Dt(garch_rv_sb[[1]]$est_vol_in_s,
                           garch_rv_sb[[2]]$est_vol_in_s, db)
ed_table11 <- build_eps_Dt(garch_rvx_no_sb[[1]]$est_vol_in_s,
                           garch_rvx_no_sb[[2]]$est_vol_in_s, db)
ed_table12 <- build_eps_Dt(garch_rvx_sb[[1]]$est_vol_in_s,
                           garch_rvx_sb[[2]]$est_vol_in_s, db)

# ---
dcc9 <- estimate_dcc_gen(ll_fun = dccmidas_ll_rc_no_sb,
                          mat_fun = dccmidas_mat_est_rc_no_sb,
                          eps_t = ed_table9$eps_t, D_t = ed_table9$D_t, TT = TT,
                          N_c = N_c, K_c = K_corr,
                          has_X = FALSE, dummy_dcc = NULL)

# --- 
dcc10 <- estimate_dcc_gen(ll_fun = dccmidas_ll_rc_sb,
                           mat_fun = dccmidas_mat_est_rc_sb,
                           eps_t = ed_table10$eps_t, D_t = ed_table10$D_t, TT = TT,
                           N_c = N_c, K_c = K_corr,
                           has_X = FALSE, dummy_dcc = dummy_dcc_input,
                           begin_val_override = coef(dcc9$m_est))

# ---
dcc11 <- estimate_dcc_gen(ll_fun = dccmidas_ll,
                           mat_fun = dccmidas_mat_est,
                           eps_t = ed_table11$eps_t, D_t = ed_table11$D_t, TT = TT,
                           N_c = N_c, K_c = K_corr, dcc_mv = mv_x,
                           has_X = TRUE, dummy_dcc = NULL)

# --- 
dcc12 <- estimate_dcc_gen(ll_fun = dccmidas_ll_sb,
                           mat_fun = dccmidas_mat_est_sb,
                           eps_t = ed_table12$eps_t, D_t = ed_table12$D_t, TT = TT,
                           N_c = N_c, K_c = K_corr, dcc_mv = mv_x,
                           has_X = TRUE, dummy_dcc = dummy_dcc_input,
                           begin_val_override = coef(dcc11$m_est))

cat("\n\n")
cat(strrep("#", 80), "\n")
cat("  ALL 12 RESULT TABLES: BTC -", toupper(INDEX_TICKER), "\n")
cat(strrep("#", 80), "\n")

# ---- GARCH-MIDAS: RV only ----
cat("\n", strrep("=", 60), "\n")
cat(" [Table 1] GARCH-MIDAS (RV only, no SB) — BTC\n")
cat(strrep("=", 60), "\n")
print(garch_rv_no_sb[[1]]$rob_coef_mat)
cat("LogL:", garch_rv_no_sb[[1]]$loglik,
    " AIC:", garch_rv_no_sb[[1]]$inf_criteria[1],
    " BIC:", garch_rv_no_sb[[1]]$inf_criteria[2], "\n")

cat("\n", strrep("=", 60), "\n")
cat(" [Table 2] GARCH-MIDAS (RV only, no SB) —", toupper(INDEX_TICKER), "\n")
cat(strrep("=", 60), "\n")
print(garch_rv_no_sb[[2]]$rob_coef_mat)
cat("LogL:", garch_rv_no_sb[[2]]$loglik,
    " AIC:", garch_rv_no_sb[[2]]$inf_criteria[1],
    " BIC:", garch_rv_no_sb[[2]]$inf_criteria[2], "\n")

cat("\n", strrep("=", 60), "\n")
cat(" [Table 3] GARCH-MIDAS (RV only, with SB) — BTC\n")
cat(strrep("=", 60), "\n")
print(garch_rv_sb[[1]]$rob_coef_mat)
cat("LogL:", garch_rv_sb[[1]]$loglik,
    " AIC:", garch_rv_sb[[1]]$inf_criteria[1],
    " BIC:", garch_rv_sb[[1]]$inf_criteria[2], "\n")

cat("\n", strrep("=", 60), "\n")
cat(" [Table 4] GARCH-MIDAS (RV only, with SB) —", toupper(INDEX_TICKER), "\n")
cat(strrep("=", 60), "\n")
print(garch_rv_sb[[2]]$rob_coef_mat)
cat("LogL:", garch_rv_sb[[2]]$loglik,
    " AIC:", garch_rv_sb[[2]]$inf_criteria[1],
    " BIC:", garch_rv_sb[[2]]$inf_criteria[2], "\n")

# ---- GARCH-MIDAS-X: RV + GEPU ----
cat("\n", strrep("=", 60), "\n")
cat(" [Table 5] GARCH-MIDAS-X (RV+GEPU, no SB) — BTC\n")
cat(strrep("=", 60), "\n")
print(garch_rvx_no_sb[[1]]$rob_coef_mat)
cat("LogL:", garch_rvx_no_sb[[1]]$loglik,
    " AIC:", garch_rvx_no_sb[[1]]$inf_criteria[1],
    " BIC:", garch_rvx_no_sb[[1]]$inf_criteria[2], "\n")

cat("\n", strrep("=", 60), "\n")
cat(" [Table 6] GARCH-MIDAS-X (RV+GEPU, no SB) —", toupper(INDEX_TICKER), "\n")
cat(strrep("=", 60), "\n")
print(garch_rvx_no_sb[[2]]$rob_coef_mat)
cat("LogL:", garch_rvx_no_sb[[2]]$loglik,
    " AIC:", garch_rvx_no_sb[[2]]$inf_criteria[1],
    " BIC:", garch_rvx_no_sb[[2]]$inf_criteria[2], "\n")

cat("\n", strrep("=", 60), "\n")
cat(" [Table 7] GARCH-MIDAS-X (RV+GEPU, with SB) — BTC\n")
cat(strrep("=", 60), "\n")
print(garch_rvx_sb[[1]]$rob_coef_mat)
cat("LogL:", garch_rvx_sb[[1]]$loglik,
    " AIC:", garch_rvx_sb[[1]]$inf_criteria["AIC"],
    " BIC:", garch_rvx_sb[[1]]$inf_criteria["BIC"], "\n")

cat("\n", strrep("=", 60), "\n")
cat(" [Table 8] GARCH-MIDAS-X (RV+GEPU, with SB) —", toupper(INDEX_TICKER), "\n")
cat(strrep("=", 60), "\n")
print(garch_rvx_sb[[2]]$rob_coef_mat)
cat("LogL:", garch_rvx_sb[[2]]$loglik,
    " AIC:", garch_rvx_sb[[2]]$inf_criteria["AIC"],
    " BIC:", garch_rvx_sb[[2]]$inf_criteria["BIC"], "\n")

# ---- DCC tables ----
for (k in 9:12) {
  obj <- get(paste0("dcc", k))
  cat("\n", strrep("=", 60), "\n")
  cat(sprintf(" [Table %d] DCC-MIDAS%s%s\n",
              k,
              if (obj$has_X) "-X" else "",
              if (obj$n_breaks > 0) ", with SB" else ", no SB"))
  cat(strrep("=", 60), "\n")
  print(obj$coef_mat)
  cat("LogL:", obj$llk, " AIC:", obj$AIC, " BIC:", obj$BIC,
      "  n_params:", obj$n_params, "  n_breaks:", obj$n_breaks, "\n")
}

# 4. HEDGING EFFECTIVENESS  ===========================================================================================================================                                                         
db_no_xts <- zoo::coredata(db)
R_B_raw   <- db_no_xts[, 1]
R_I_raw   <- db_no_xts[, 2]
var_I     <- var(R_I_raw, na.rm = TRUE)

compute_HE <- function(H_t, R_B, R_I, var_I) {
  TT    <- dim(H_t)[3]
  gamma <- sapply(1:TT, function(tt) H_t[1, 2, tt] / H_t[1, 1, tt])
  R_H   <- R_I - gamma * R_B
  var_H <- var(R_H, na.rm = TRUE)
  list(gamma      = gamma,
       mean_gamma = mean(gamma, na.rm = TRUE),
       sd_gamma   = sd(gamma, na.rm = TRUE),
       R_H        = R_H,
       var_H      = var_H,
       HE         = (var_I - var_H) / var_I)
}

he9  <- compute_HE(dcc9$matrices$H_t,  R_B_raw, R_I_raw, var_I)
he10 <- compute_HE(dcc10$matrices$H_t, R_B_raw, R_I_raw, var_I)
he11 <- compute_HE(dcc11$matrices$H_t, R_B_raw, R_I_raw, var_I)
he12 <- compute_HE(dcc12$matrices$H_t, R_B_raw, R_I_raw, var_I)

cat("\n=== HEDGING EFFECTIVENESS COMPARISON ===\n")
cat(sprintf("%-35s %8s %8s %10s\n", "Model", "Mean OHR", "SD OHR", "HE"))
cat(strrep("-", 65), "\n")
cat(sprintf("%-35s %8.4f %8.4f %10.6f\n", "Table 9  DCC-RC  (no SB)",
            he9$mean_gamma,  he9$sd_gamma,  he9$HE))
cat(sprintf("%-35s %8.4f %8.4f %10.6f\n", "Table 10 DCC-RC  (with SB)",
            he10$mean_gamma, he10$sd_gamma, he10$HE))
cat(sprintf("%-35s %8.4f %8.4f %10.6f\n", "Table 11 DCC-RCX (no SB)",
            he11$mean_gamma, he11$sd_gamma, he11$HE))
cat(sprintf("%-35s %8.4f %8.4f %10.6f\n", "Table 12 DCC-RCX (with SB)",
            he12$mean_gamma, he12$sd_gamma, he12$HE))
cat(strrep("-", 65), "\n")
cat(sprintf("Delta HE (GEPU effect, no SB):   %+.6f\n", he11$HE - he9$HE))
cat(sprintf("Delta HE (SB effect, no GEPU):   %+.6f\n", he10$HE - he9$HE))
cat(sprintf("Delta HE (full extension):       %+.6f\n", he12$HE - he9$HE))

cat("\n=== GARCH-MIDAS LR TESTS (H0: theta_X = 0) ===\n")
lr_test <- function(ll_unrest, ll_rest, df, label) {
  stat <- 2 * (ll_unrest - ll_rest)
  pval <- pchisq(stat, df = df, lower.tail = FALSE)
  cat(sprintf("%-55s LR=%8.4f  df=%d  p=%.6f %s\n", label, stat, df, pval,
              if (pval < 0.01) "***" else if (pval < 0.05) "**" else if (pval < 0.10) "*" else ""))
}

lr_test(garch_rvx_no_sb[[1]]$loglik, garch_rv_no_sb[[1]]$loglik, 2,
        "BTC:  GARCH-MIDAS-X vs GARCH-MIDAS (no SB)")
lr_test(garch_rvx_no_sb[[2]]$loglik, garch_rv_no_sb[[2]]$loglik, 2,
        paste0(toupper(INDEX_TICKER), ": GARCH-MIDAS-X vs GARCH-MIDAS (no SB)"))

lr_test(garch_rvx_sb[[1]]$loglik, garch_rvx_no_sb[[1]]$loglik,
        garch_rvx_sb[[1]]$n_params - garch_rvx_no_sb[[1]]$n_params,
        "BTC:  GARCH-MIDAS-X+SB vs GARCH-MIDAS-X (SB added)")
lr_test(garch_rvx_sb[[2]]$loglik, garch_rvx_no_sb[[2]]$loglik,
        garch_rvx_sb[[2]]$n_params - garch_rvx_no_sb[[2]]$n_params,
        paste0(toupper(INDEX_TICKER), ": GARCH-MIDAS-X+SB vs GARCH-MIDAS-X (SB added)"))

cat("\n=== WALD JOINT TESTS (H0: all break coefficients = 0) ===\n")
wald_joint_test <- function(garch_obj, break_names, label) {
  coef_mat <- garch_obj$rob_coef_mat
  idx <- match(break_names, rownames(coef_mat))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) { cat(sprintf("  %-50s  No breaks found\n", label)); return() }
  theta_hat <- coef_mat$Estimate[idx]
  se_hat    <- coef_mat$`Std. Error`[idx]
  valid     <- !is.na(se_hat) & se_hat > 0
  wald_stat <- sum((theta_hat[valid] / se_hat[valid])^2)
  df_wald   <- sum(valid)
  p_wald    <- pchisq(wald_stat, df = df_wald, lower.tail = FALSE)
  cat(sprintf("  %-50s W=%8.4f  df=%d  p=%.6f %s\n", label, wald_stat, df_wald, p_wald,
              if (p_wald < 0.01) "***" else if (p_wald < 0.05) "**" else if (p_wald < 0.10) "*" else ""))
}
rv_brk <- paste0("theta_RV_break_", 1:n_breaks)
x_brk  <- paste0("theta_X_break_",  1:n_breaks)
wald_joint_test(garch_rv_sb[[1]],  rv_brk,          "BTC  Table 3: RV breaks")
wald_joint_test(garch_rv_sb[[2]],  rv_brk,          paste0(toupper(INDEX_TICKER), " Table 4: RV breaks"))
wald_joint_test(garch_rvx_sb[[1]], rv_brk,          "BTC  Table 7: RV breaks (with GEPU)")
wald_joint_test(garch_rvx_sb[[1]], x_brk,           "BTC  Table 7: GEPU breaks")
wald_joint_test(garch_rvx_sb[[1]], c(rv_brk, x_brk),"BTC  Table 7: ALL breaks (joint)")
wald_joint_test(garch_rvx_sb[[2]], rv_brk,          paste0(toupper(INDEX_TICKER), " Table 8: RV breaks (with GEPU)"))
wald_joint_test(garch_rvx_sb[[2]], x_brk,           paste0(toupper(INDEX_TICKER), " Table 8: GEPU breaks"))
wald_joint_test(garch_rvx_sb[[2]], c(rv_brk, x_brk),paste0(toupper(INDEX_TICKER), " Table 8: ALL breaks (joint)"))

cat("\n=== DCC-MIDAS LR TESTS ===\n")
cat("NOTE: Each DCC model uses matched first-stage GARCH residuals.\n")
cat("      LR tests are informative but not strictly nested due to different first-stages.\n\n")
lr_test(dcc11$llk, dcc9$llk,  dcc11$n_params - dcc9$n_params,  "Table11 vs Table9  (GEPU effect, no SB)")
lr_test(dcc10$llk, dcc9$llk,  dcc10$n_params - dcc9$n_params,  "Table10 vs Table9  (SB effect, no GEPU)")
lr_test(dcc12$llk, dcc11$llk, dcc12$n_params - dcc11$n_params, "Table12 vs Table11 (SB added to GEPU)")

cat("\n=== FORECAST COMPARISON TESTS ===\n")
loss_base <- he9$R_H^2
loss_full <- he12$R_H^2
valid <- !is.na(loss_base) & !is.na(loss_full)
dm_result <- tryCatch(
  forecast::dm.test(ts(loss_base[valid]), ts(loss_full[valid]),
                    alternative = "greater", h = 1, power = 1),
  error = function(e) { cat("DM test error:", e$message, "\n"); NULL })
if (!is.null(dm_result)) {
  cat(sprintf("DM Test  (Table12 vs Table9): stat=%7.4f  p=%.4f %s\n",
              dm_result$statistic, dm_result$p.value,
              if (dm_result$p.value < 0.05) "**" else ""))
}

cw_diff  <- (loss_base - loss_full) + (he9$R_H - he12$R_H)^2
valid_cw <- !is.na(cw_diff)
cw_mean  <- mean(cw_diff[valid_cw])
cw_nw_se <- tryCatch({
  
  T_cw <- sum(valid_cw)
  bw   <- floor(T_cw^(1/3))
  cw_centered <- cw_diff[valid_cw] - cw_mean
  gamma0 <- mean(cw_centered^2)
  gammaJ <- sapply(1:bw, function(j) mean(cw_centered[1:(T_cw-j)] * cw_centered[(j+1):T_cw]))
  nw_var <- gamma0 + 2 * sum((1 - (1:bw)/(bw+1)) * gammaJ)
  sqrt(nw_var / T_cw)
}, error = function(e) sd(cw_diff[valid_cw]) / sqrt(sum(valid_cw)))
cw_stat  <- cw_mean / cw_nw_se
cw_pval  <- 1 - pnorm(cw_stat)
cat(sprintf("CW Test  (Table12 vs Table9): stat=%7.4f  p=%.4f %s\n",
            cw_stat, cw_pval,
            if (cw_pval < 0.05) "**" else ""))

loss_gepu <- he11$R_H^2
cw_diff2  <- (loss_base - loss_gepu) + (he9$R_H - he11$R_H)^2
valid_cw2 <- !is.na(cw_diff2)
cw_mean2  <- mean(cw_diff2[valid_cw2])
cw_se2    <- sd(cw_diff2[valid_cw2]) / sqrt(sum(valid_cw2))
cw_stat2  <- cw_mean2 / cw_se2
cw_pval2  <- 1 - pnorm(cw_stat2)
cat(sprintf("CW Test  (Table11 vs Table9): stat=%7.4f  p=%.4f %s\n",
            cw_stat2, cw_pval2,
            if (cw_pval2 < 0.05) "**" else ""))

cat("\n=== RESIDUAL DIAGNOSTICS ===\n")

diag_one <- function(garch_obj, r_t_i, label) {
  vol_total <- as.numeric(garch_obj$est_vol_in_s)
  resid_raw <- as.numeric(r_t_i)

  mu_hat <- garch_obj$rob_coef_mat["mu", "Estimate"]
  xi <- (resid_raw - mu_hat) / vol_total
  xi_clean <- xi[!is.na(xi) & is.finite(xi)]

  lb10   <- Box.test(xi_clean,   lag = 10, type = "Ljung-Box")
  lb20   <- Box.test(xi_clean,   lag = 20, type = "Ljung-Box")
  lb10sq <- Box.test(xi_clean^2, lag = 10, type = "Ljung-Box")
  lb20sq <- Box.test(xi_clean^2, lag = 20, type = "Ljung-Box")

  arch_lm <- tryCatch({
    xi2 <- xi_clean^2
    emb <- embed(xi2, 11)
    lm_fit <- lm(emb[,1] ~ emb[,-1])
    r2 <- summary(lm_fit)$r.squared
    n_lm <- nrow(emb)
    stat_lm <- n_lm * r2
    p_lm <- pchisq(stat_lm, df = 10, lower.tail = FALSE)
    list(stat = stat_lm, p = p_lm)
  }, error = function(e) list(stat = NA, p = NA))

  cat(sprintf("  %-40s Q10=%.2f(%.3f) Q20=%.2f(%.3f) Q10sq=%.2f(%.3f) Q20sq=%.2f(%.3f) LM=%.2f(%.3f)\n",
              label,
              lb10$statistic, lb10$p.value,
              lb20$statistic, lb20$p.value,
              lb10sq$statistic, lb10sq$p.value,
              lb20sq$statistic, lb20sq$p.value,
              ifelse(is.na(arch_lm$stat), NA, arch_lm$stat),
              ifelse(is.na(arch_lm$p), NA, arch_lm$p)))
}

cat("\n  BTC:\n")
diag_one(garch_rv_no_sb[[1]],  r_t_[[1]], "Table 1: RV only, no SB")
diag_one(garch_rv_sb[[1]],     r_t_[[1]], "Table 3: RV only, with SB")
diag_one(garch_rvx_no_sb[[1]], r_t_[[1]], "Table 5: RV+GEPU, no SB")
diag_one(garch_rvx_sb[[1]],    r_t_[[1]], "Table 7: RV+GEPU, with SB (full)")

cat(sprintf("\n  %s:\n", toupper(INDEX_TICKER)))
diag_one(garch_rv_no_sb[[2]],  r_t_[[2]], "Table 2: RV only, no SB")
diag_one(garch_rv_sb[[2]],     r_t_[[2]], "Table 4: RV only, with SB")
diag_one(garch_rvx_no_sb[[2]], r_t_[[2]], "Table 6: RV+GEPU, no SB")
diag_one(garch_rvx_sb[[2]],    r_t_[[2]], "Table 8: RV+GEPU, with SB (full)")

# TRANSACTION COST ANALYSIS                                                
cat("\n=== HEDGING EFFECTIVENESS WITH TRANSACTION COSTS ===\n")
compute_HE_tc <- function(H_t, R_B, R_I, var_I, cost_bps) {
  c <- cost_bps / 10000 
  TT <- dim(H_t)[3]
  gamma <- sapply(1:TT, function(tt) H_t[1, 2, tt] / H_t[1, 1, tt])
  turnover <- c(0, abs(diff(gamma)))
  R_H_net <- R_I - gamma * R_B - c * turnover * abs(R_B)
  var_H_net <- var(R_H_net, na.rm = TRUE)
  (var_I - var_H_net) / var_I
}

tc_levels <- c(0, 10, 25, 50)
cat(sprintf("%-35s %8s %8s %8s %8s\n", "Model", "c=0bp", "c=10bp", "c=25bp", "c=50bp"))
cat(strrep("-", 75), "\n")
for (k in c(9, 10, 11, 12)) {
  H_t_k <- get(paste0("dcc", k))$matrices$H_t
  label <- paste0("Table ", k, " ",
                  if (k %in% c(11,12)) "DCC-X" else "DCC-RC",
                  if (k %in% c(10,12)) "+SB" else "")
  he_vals <- sapply(tc_levels, function(tc)
    compute_HE_tc(H_t_k, R_B_raw, R_I_raw, var_I, tc))
  cat(sprintf("%-35s %8.6f %8.6f %8.6f %8.6f\n", label,
              he_vals[1], he_vals[2], he_vals[3], he_vals[4]))
}

cat("\n\n========== ANALYSIS COMPLETE ==========\n")
cat(sprintf("Pair: BTC - %s | Sample: %s to %s\n",
            toupper(INDEX_TICKER), start_date_is, end_date_is))
cat(sprintf("K_vol=%d, K_corr=%d, N_c=%d, n_breaks=%d\n", K_vol, K_corr, N_c, n_breaks))