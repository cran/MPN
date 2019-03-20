# Point estimate for MPN -------------------------------------------------------

.ptEst_MPN <- function(positive, tubes, amount) {
  # https://www.fda.gov/Food/FoodScienceResearch/LaboratoryMethods/ucm109656.htm
  .scoreFnc <- function(lambda, positive, tubes, amount) {
    #score function: Find root to maximize likelihood.
    LHS <- sum((positive * amount) / (1 - exp(-lambda * amount)))
    RHS <- sum(tubes * amount)
    LHS - RHS
  }

  all_negative <- sum(positive) == 0
  all_positive <- identical(positive, tubes)

  if (all_negative) {
    return(0)
  } else if (all_positive) {
    return(Inf)
  } else {
    MPN <- uniroot(.scoreFnc, interval = c(1e-04, 1e+04),
                   positive = positive, tubes = tubes, amount = amount,
                   extendInt = "downX", maxiter = 1e+04)$root
    return(MPN)
  }
}

.ptEstAdj_MPN <- function(MPN, tubes, amount) {
  #bias adjustment for MPN (Salama et al., 1978; Haas, 1989)
  if (MPN == 0) {
    return(0)
  } else if (is.infinite(MPN)) {
    return(NA)
  } else {
    amount2 <- amount ^ 2
    amount3 <- amount ^ 3
    lambda_V <- MPN * amount
    one_or_more <- 1 - exp(-lambda_V)
    zi  <- tubes * one_or_more
    cosh_term <- cosh(lambda_V) - 1
    D   <- sum((amount2 * zi) / (2 * cosh_term))
    wi1 <- amount2 / (2 * (one_or_more ^ 2) * (D ^ 3))
    wi2 <- (amount3 * zi * sinh(lambda_V)) / (cosh_term ^ 2)
    wi2[is.nan(wi2)] <- 0  #Inf/Inf is NaN (denom goes to infinity more quickly)
    wi2 <- sum(wi2)
    wi3 <- amount3 / (one_or_more * cosh_term * (D ^ 2))
    wi  <- wi1 * wi2 - wi3
    return(MPN - 0.5 * sum(wi * zi * exp(-lambda_V)))
  }
}

# Likelihood, LRT, etc. --------------------------------------------------------

# See Ridout (1994) - "A Comparison of CI Methods..."
# Section 2.2 "Intervals Based on the LRT"

.logL_MPN <- function(lambda, positive, tubes, amount) {
  #log-likelihood
  lambda_amount <- lambda * amount
  first_term    <- positive * log(1 - exp(-lambda_amount))
  second_term   <- (tubes - positive) * lambda_amount
  sum(first_term - second_term)
}

.logLR_MPN <- function(lambda, lambda_hat, positive, tubes, amount) {
  #log-likelihood ratio
  logL_alt  <- .logL_MPN(lambda_hat, positive, tubes, amount)
  logL_null <- .logL_MPN(lambda, positive, tubes, amount)
  2 * (logL_alt - logL_null)
}

.logLRroot_MPN <- function(lambda, lambda_hat, positive, tubes, amount,
                           crit_val) {
  #Find roots to get LR confidence limits
  log_like_ratio <- .logLR_MPN(lambda, lambda_hat, positive, tubes, amount)
  log_like_ratio - crit_val
}


# Confidence intervals ---------------------------------------------------------

.jarvisCI_MPN <- function(MPN, positive, tubes, amount, conf_level) {
  # See Jarvis (2010) - "Reconsideration of the derivation of..."
  all_negative <- sum(positive) == 0
  all_positive <- identical(positive, tubes)
  var_logMPN <- NA
  variance   <- NA
  sig_level  <- 1 - conf_level

  if (!all_negative && !all_positive) {
    exp_term <- exp(-amount * MPN)
    numer    <- positive * amount ^ 2 * exp_term
    denom    <- (1 - exp_term) ^ 2
    variance <- 1 / sum(numer / denom)
  }

  if (all_negative) {
    LB <- 0
    UB <- log(1 / sig_level) / sum(amount * tubes)
  } else if (all_positive) {
    .fLB <- function(LB) {
      log(1 / sig_level) + sum(tubes * log(1 - exp(-amount * LB)))
    }
    LB <- uniroot(.fLB, interval = c(1e-02, 1e+03),
                  extendInt = "upX", maxiter = 1e+04)$root
    UB <- Inf
  } else {
    crit_val   <- qnorm(sig_level / 2, lower.tail = FALSE)  #asym. normal
    var_logMPN <- variance / (MPN ^ 2)  #delta method
    SE_log     <- sqrt(var_logMPN)
    ME_log     <- crit_val * SE_log
    LB <- MPN * exp(-ME_log)
    UB <- MPN * exp(ME_log)
  }

  list(variance = variance, var_logMPN = var_logMPN, LB = LB, UB = UB)
}

.likeRatioCI_MPN <- function(MPN, positive, tubes, amount, conf_level) {
  #Likelihood ratio confidence limits
  #See Ridout (1994) - "A Comparison of CI Methods..."
  all_negative <- sum(positive) == 0
  all_positive <- identical(positive, tubes)
  crit_val <- qchisq(conf_level, df = 1, lower.tail = TRUE)
  if (all_negative || all_positive) {
    jarvis_bounds <- .jarvisCI_MPN(MPN, positive, tubes, amount, conf_level)
    LB <- jarvis_bounds$LB
    UB <- jarvis_bounds$UB
  } else {
    LB <- uniroot(.logLRroot_MPN, interval = c(1e-04, MPN),
                  lambda_hat = MPN, positive = positive, tubes = tubes,
                  amount = amount, crit_val = crit_val,
                  extendInt  = "no", maxiter = 1e+04)$root
    UB <- uniroot(.logLRroot_MPN, interval = c(MPN, 5 * MPN),
                  lambda_hat = MPN, positive = positive, tubes = tubes,
                  amount = amount, crit_val = crit_val,
                  extendInt  = "upX", maxiter = 1e+04)$root
  }
  list(LB = LB, UB = UB)
}


# Blodgett's rarity index ------------------------------------------------------

.rarity_MPN <- function(MPN, positive, tubes, amount) {

  all_negative <- sum(positive) == 0
  all_positive <- identical(positive, tubes)

  if (all_negative || all_positive) {
    return(1)
  } else {
    probs        <- 1 - exp(-MPN * amount)
    positive_max <- floor(probs * (tubes + 1))  #Jarvis p.1664
    positive_max <- pmin(positive_max, tubes)   # and Blodgett (2005) App.C
    most_prob_L  <- dbinom(x = positive_max, size = tubes, prob = probs)
    actual_L     <- dbinom(x = positive, size = tubes, prob = probs)
    rarity       <- prod(actual_L) / prod(most_prob_L)
    return(rarity)
  }
}
