# Point estimate for APC -------------------------------------------------------

.ptEst_APC <- function(count, amount_scor, amount_tntc, tntc_limit) {
  if (sum(count) == 0) {
    return(0)
  }
  lower_search <- 1e-5
  upper_search <- 2 * tntc_limit / min(1, amount_scor, amount_tntc)
  if (length(amount_tntc) > 0) {
    lower_search <- optimize(.logLscor_APC, count = count,
                             amount_scor = amount_scor,
                             lower = lower_search, upper = upper_search,
                             tol = 1e-05, maximum = TRUE)$maximum
  }
  optimize(.logL_APC, count = count, amount_scor = amount_scor,
           amount_tntc = amount_tntc, tntc_limit = tntc_limit,
           lower = lower_search, upper = upper_search, tol = 1e-05,
           maximum = TRUE)$maximum
}

# Likelihood, LRT, etc. --------------------------------------------------------

# See Haas (2014) - "Quantitative Microbial Risk Assessment..." (Ch.6)

.logLscor_APC <- function(lambda, count, amount_scor) {
  # The scorable portion of the log-likelihood.
  # Eq 6.8 (p.162)
  scorable1 <- sum(count * log(lambda * amount_scor))
  scorable2 <- lambda * sum(amount_scor)
  scorable1 - scorable2
}

.logL_APC <- function(lambda, count, amount_scor, amount_tntc, tntc_limit) {
  # Entire likelihood (scorable & TNTC).
  # Eq 6.11 (p.163)
  # See Eq. 6.6 (p.161) & help(pgamma) for the incomplete gamma relationship.
  scorable <- .logLscor_APC(lambda, count = count, amount_scor = amount_scor)
  if (length(amount_tntc) > 0) {
    incomplete_gamma <- pgamma(tntc_limit - 1, lambda * amount_tntc)
    tntc <- sum(log(1 - incomplete_gamma))
    if (is.infinite(tntc)) {
      error_msg <- paste("one or more of these TNTC plates are extremely",
                         "unlikely\nconsidering the observed scorable plates")
      stop(error_msg)
    }
    return(scorable + tntc)
  } else {
    return(scorable)
  }
}

.logLR_APC <- function(lambda, lambda_hat, count, amount_scor, amount_tntc,
                       tntc_limit) {
  # See Eq.6.38 (p.185)
  scorable1 <- count * log(lambda / lambda_hat)
  scorable2 <- amount_scor * (lambda - lambda_hat)
  scorable  <- sum(scorable1 - scorable2)
  if (length(amount_tntc) > 0) {
    tntc_sum_limit <- tntc_limit - 1
    incomplete_gamma1 <- 1 - pgamma(tntc_sum_limit, lambda * amount_tntc)
    incomplete_gamma2 <- 1 - pgamma(tntc_sum_limit, lambda_hat * amount_tntc)
    tntc <- sum(log(incomplete_gamma1 / incomplete_gamma2))
    return(-2 * (scorable + tntc))
  } else {
    return(-2 * scorable)
  }
}

.logLRroot_APC <- function(lambda, lambda_hat, count, amount_scor, amount_tntc,
                           tntc_limit, crit_val) {
  # Find roots to identify LRT CI limits.
  like_ratio <- .logLR_APC(lambda, lambda_hat = lambda_hat, count = count,
                           amount_scor = amount_scor, amount_tntc = amount_tntc,
                           tntc_limit = tntc_limit)
  like_ratio - crit_val
}

# Confidence intervals ---------------------------------------------------------

.likeRatioCI_APC <- function(lambda_hat, count, amount_scor, amount_tntc,
                             tntc_limit, conf_level) {
  if (sum(count) == 0) {
    LB <- UB <- NA
  } else {
    crit_val <- qchisq(conf_level, df = 1, lower.tail = TRUE)
    LB <- uniroot(.logLRroot_APC, interval = c(.5 * lambda_hat, lambda_hat),
                  lambda_hat = lambda_hat, count = count,
                  amount_scor = amount_scor, amount_tntc = amount_tntc,
                  tntc_limit = tntc_limit, crit_val = crit_val,
                  extendInt = "downX", maxiter = 1e+04)$root

    UB <- uniroot(.logLRroot_APC, interval = c(lambda_hat, 2 * lambda_hat),
                  lambda_hat = lambda_hat, count = count,
                  amount_scor = amount_scor, amount_tntc = amount_tntc,
                  tntc_limit = tntc_limit, crit_val = crit_val,
                  extendInt = "upX", maxiter = 1e+04)$root
  }
  list(LB = LB, UB = UB)
}
