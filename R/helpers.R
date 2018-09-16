.ptEstimate <- function(positive, tubes, amount) {

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
                   extendInt = "yes", maxiter = 1e+04)$root
    return(MPN)
  }
}


.jarvisCI <- function(MPN, positive, tubes, amount, sig_level) {

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
    LB <- uniroot(.fLB, c(0, 999999))$root
    UB <- Inf
  } else {
    crit_val   <- qnorm(sig_level / 2, lower.tail = FALSE)  #asym. normal
    var_logMPN <- variance / (MPN ^ 2)  #delta method
    SE_log     <- sqrt(var_logMPN)
    ME_log     <- crit_val * SE_log
    LB <- MPN * exp(-ME_log)
    UB <- MPN * exp(ME_log)
  }

  list(variance = variance, var_logMPN = var_logMPN, sig_level = sig_level,
       LB = LB, UB = UB)
}


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
