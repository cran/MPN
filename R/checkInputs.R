# Purpose: Check the inputs of the functions (defensive programming)

.isMissing <- function(x) {
  #Determines if any elements in a numeric vector are missing.
  any(is.na(x) | is.nan(x) | is.infinite(x))
}

.checkInputs_mpn <- function(positive, tubes, amount, conf_level) {

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
  if (any(positive > tubes)) {
    stop("more positive tubes than possible")
  }
}
