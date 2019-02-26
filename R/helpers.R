# Point estimate for MPN -------------------------------------------------------

.ptEstimate <- function(positive, tubes, amount) {
  # https://www.fda.gov/Food/FoodScienceResearch/LaboratoryMethods/ucm109656.htm
  .scoreFnc <- function(lambda, positive, tubes, amount) {
    # Score function. Find root to maximize likelihood.
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


# Likelihood, LRT, etc. --------------------------------------------------------

# See Ridout (1994) - "A Comparison of CI Methods..."
# Section 2.2 "Intervals Based on the LRT"

.logL <- function(lambda, positive, tubes, amount) {
  #log-likelihood
  lambda_amount <- lambda * amount
  first_term    <- positive * log(1 - exp(-lambda_amount))
  second_term   <- (tubes - positive) * lambda_amount
  sum(first_term - second_term)
}

.logLR <- function(lambda, lambda_hat, positive, tubes, amount) {
  #log-likelihood ratio
  logL_alt  <- .logL(lambda_hat, positive, tubes, amount)
  logL_null <- .logL(lambda, positive, tubes, amount)
  2 * (logL_alt - logL_null)
}

.logLRroot <- function(lambda, lambda_hat, positive, tubes, amount, alpha = .05) {
  #Find roots to get LR confidence limits
  log_like_ratio <- .logLR(lambda, lambda_hat, positive, tubes, amount)
  crit_val       <- qchisq(alpha, df = 1, lower.tail = FALSE)
  log_like_ratio - crit_val
}


# Confidence intervals ---------------------------------------------------------

.jarvisCI <- function(MPN, positive, tubes, amount, sig_level) {
  # See Jarvis (2010) - "Reconsideration of the derivation of..."
  all_negative <- sum(positive) == 0
  all_positive <- identical(positive, tubes)
  var_logMPN <- NA
  variance   <- NA

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

.likeRatioCI <- function(MPN, positive, tubes, amount, sig_level) {
  #Likelihood ratio confidence limits
  #See Ridout (1994) - "A Comparison of CI Methods..."
  all_negative <- sum(positive) == 0
  all_positive <- identical(positive, tubes)
  if (all_negative || all_positive) {
    jarvis_bounds <- .jarvisCI(MPN, positive, tubes, amount, sig_level)
    LB <- jarvis_bounds$LB
    UB <- jarvis_bounds$UB
  } else {
    LB <- uniroot(.logLRroot, interval = c(1e-04, MPN),
                  lambda_hat = MPN, positive = positive, tubes = tubes,
                  amount = amount, alpha = sig_level,
                  extendInt  = "no", maxiter = 1e+04)$root
    UB <- uniroot(.logLRroot, interval = c(MPN, 5 * MPN),
                  lambda_hat = MPN, positive = positive, tubes = tubes,
                  amount = amount, alpha = sig_level,
                  extendInt  = "upX", maxiter = 1e+04)$root
  }
  list(LB = LB, UB = UB)
}


# Blodgett's rarity index ------------------------------------------------------

.rarity <- function(MPN, positive, tubes, amount) {

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
