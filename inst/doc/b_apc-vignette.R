## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----Define likelihood R function---------------------------------------------
L <- function(lambda, count, amount_scor, amount_tntc = NULL, tntc_limit = 100) {
  #likelihood
  scorable <- prod(dpois(count, lambda = lambda * amount_scor))
  if (length(amount_tntc) > 0) {
    incomplete_gamma <- pgamma(tntc_limit - 1, lambda * amount_tntc)
    tntc <- prod(1 - incomplete_gamma)
    return(scorable * tntc)
  } else {
    return(scorable)
  }
}
L_vec <- Vectorize(L, "lambda")

## ----Calculate APC estimate---------------------------------------------------
#APC calculation
library(MPN)
my_count       <- c(28, 20)          #Ni
my_amount_scor <- 1 * c(.001, .001)  #Vi
my_amount_tntc <- 1 * c(.01, .01)    #Vj
my_tntc_limit  <- c(300, 250)        #NLj
(my_apc <- apc(my_count, my_amount_scor, my_amount_tntc, my_tntc_limit))

## ----Plot likelihood and APC estimate-----------------------------------------
my_apc$APC
my_lambda <- seq(29000, 31500, length = 1000)
my_L <- L_vec(my_lambda, my_count, my_amount_scor, my_amount_tntc,
              my_tntc_limit)
plot(my_lambda, my_L, type = "l",
     ylab = "Likelihood", main = "Maximum Likelihood")
abline(v = my_apc$APC, lty = 2, col = "red")

## ----All zero counts----------------------------------------------------------
all_zero <- c(0, 0)  #Ni
(apc_all_zero <- apc(all_zero, my_amount_scor)$APC)
my_lambda <- seq(0, 1000, length = 1000)
L_all_zero  <- L_vec(my_lambda, all_zero, my_amount_scor)
plot(my_lambda, L_all_zero, type = "l", ylab = "Likelihood",
     main = "All Zeroes")
abline(v = apc_all_zero, lty = 2, col = "red")

## ----Calculate confidence interval--------------------------------------------
my_count       <- c(28, 20)
my_amount_scor <- 1 * c(.001, .001)
my_amount_tntc <- 1 * c(.01, .01)
my_tntc_limit  <- c(300, 250)
(my_apc <- apc(my_count, my_amount_scor, my_amount_tntc, my_tntc_limit))
my_lambda <- seq(25000, 36000, length = 1000)
my_L <- L_vec(my_lambda, my_count, my_amount_scor, my_amount_tntc,
              my_tntc_limit)
plot(my_lambda, my_L, type = "l",
     ylab = "Likelihood", main = "Maximum Likelihood")
abline(v = my_apc$APC, lty = 2, col = "red")
abline(v = my_apc$LB, lty = 3, col = "blue")
abline(v = my_apc$UB, lty = 3, col = "blue")

