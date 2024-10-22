# Purpose: Check the inputs of the functions (defensive programming)

.isMissing <- function(x) {
  #Determines if any elements in a numeric vector are missing.
  any(is.na(x) | is.nan(x) | is.infinite(x))
}

.checkInputs_mpn <- function(positive, tubes, amount, conf_level, tol) {

  l_positive <- length(positive)
  l_tubes    <- length(tubes)
  l_amount   <- length(amount)
  if (l_positive != l_tubes || l_tubes != l_amount) {
    stop("'positive', 'tubes', & 'amount' must be the same length")
  }
  if (.isMissing(c(positive, tubes, amount, conf_level))) {
    stop("missing values are not allowed")
  }
  if (!is.numeric(c(positive, tubes, amount, conf_level))) {
    stop("'positive', 'tubes', 'amount', & 'conf_level' must be numeric")
  }
  if (any(tubes < 1) || any(round(tubes) != tubes)) {
    stop("'tubes' must contain positive whole numbers")
  }
  if (any(positive < 0) || any(round(positive) != positive)) {
    stop("'positive' must contain non-negative whole numbers")
  }
  if (any(amount <= 0)) {
    stop("'amount' must contain positive values")
  }
  if (length(amount) > 1 && max(diff(amount)) >= 0) {
    stop("'amount' must be in descending order")
  }
  if (length(conf_level) != 1) {
    stop("'conf_level' must have length of 1")
  }
  if (conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be between 0 & 1")
  }
  if (length(tol) != 1) {
    stop("'tol' must have length of 1")
  }
  if (tol <= 0 || tol > 1e-3) {
    stop("'tol' must be positive and at most 1e-3")
  }
  if (any(positive > tubes)) {
    stop("more positive tubes than possible")
  }
}

.checkInputs_apc <- function(count, amount_scor, amount_tntc, tntc_limit,
                             conf_level, tol) {

  l_count       <- length(count)
  l_amount_scor <- length(amount_scor)
  if (l_count != l_amount_scor) {
    stop("'count' &  'amount_scor' must be the same length")
  }
  if (l_count < 1) {
    stop("there must be at least one scorable plate")
  }
  if (any(is.na(count)) || any(is.na(amount_scor))) {
    stop("missing values not allowed for 'count' or 'amount_scor'")
  }
  if (any(count < 0) || any(round(count) != count)) {
    stop("'count' must contain non-negative whole numbers")
  }
  if (any(amount_scor <= 0)) {
    stop("'amount_scor' must contain positive values")
  }
  if (!is.null(amount_tntc) && any(amount_tntc <= 0)) {
    stop("'amount_tntc' must be NULL or contain positive values")
  }
  if (length(conf_level) != 1) {
    stop("'conf_level' must have length of 1")
  }
  if (conf_level <= 0 || conf_level >= 1) {
    stop("'conf_level' must be between 0 & 1")
  }
  if (length(tol) != 1) {
    stop("'tol' must have length of 1")
  }
  if (tol <= 0 || tol > 1e-3) {
    stop("'tol' must be positive and at most 1e-3")
  }
  if (sum(count) == 0 && length(amount_tntc) != 0) {
    stop("if all scorable counts are zero, there should not be any TNTC plates")
  }
}
