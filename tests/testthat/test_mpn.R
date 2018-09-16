context("mpn function - verify results")
#Test if the function mpn() works properly.

# Quick comparison table -------------------------------------------------------

test_that("APHA quick table - MPN point estimates correct", {
  tubes1 <- c(3, 3, 3)
  amount <- c(10, 1, .1)
  expect_equal(
    signif(mpn(positive = c(1, 1, 1), tubes = tubes1, amount = amount)$MPN, 2),
    .11
  )
  expect_equal(
    signif(mpn(positive = c(2, 1, 0), tubes = tubes1, amount = amount)$MPN, 2),
    .15
  )
  expect_equal(
    signif(mpn(positive = c(3, 1, 0), tubes = tubes1, amount = amount)$MPN, 2),
    .43
  )
  tubes2 <- c(5, 5, 5)
  expect_equal(
    signif(mpn(positive = c(1, 1, 1), tubes = tubes2, amount = amount)$MPN, 2),
    .061
  )
  expect_equal(
    signif(mpn(positive = c(2, 1, 0), tubes = tubes2, amount = amount)$MPN, 2),
    .068
  )
  expect_equal(
    signif(mpn(positive = c(3, 1, 0), tubes = tubes2, amount = amount)$MPN, 2),
    .11
  )
  expect_equal(
    signif(mpn(positive = c(5, 2, 0), tubes = tubes2, amount = amount)$MPN, 2),
    .49
  )
})


# Blodgett's articles ----------------------------------------------------------

test_that("Blodgett - MPN & RI match spreadsheet", {
  amount <- c(.1, .01, .001)
  x <- mpn(positive = c(5, 3, 0), tubes = c(10, 10, 10), amount = amount)
  expect_equal(x$MPN, 9.950764229, tol = .001)
  expect_equal(x$RI, .0924576, tol = .001)

  y <- mpn(positive = c(4, 2, 2), tubes = c(5, 5, 4), amount = amount)
  expect_equal(y$MPN, 31.918, tol = .002)
  expect_equal(y$RI, 1.02e-3, tol = .0005)

  z <- mpn(positive = c(3, 0, 2), tubes = c(7, 5, 7), amount = amount)
  expect_equal(z$MPN, 8.549, tol = .002)
  expect_equal(z$RI, 1.15e-3, tol = .0005)

  q <- mpn(positive = c(2, 3, 4), tubes = c(9, 3, 4), amount = amount)
  expect_equal(q$MPN, 11.274, tol = .002)
  expect_equal(q$RI, 6.34e-13, tol = 1e-12)
})

test_that("Blodgett - MPN matches Table 1, Blodgett (2002)", {
  tubes  <- c(5, 5)
  amount <- c(1, .1)
  expect_equal(mpn(positive = c(1, 5), tubes = tubes, amount = amount)$MPN,
               1.29, tol = .01)
  expect_equal(mpn(positive = c(2, 4), tubes = tubes, amount = amount)$MPN,
               1.48, tol = .01)
  expect_equal(mpn(positive = c(3, 3), tubes = tubes, amount = amount)$MPN,
               1.75, tol = .01)
  expect_equal(mpn(positive = c(4, 2), tubes = tubes, amount = amount)$MPN,
               2.21, tol = .01)
  expect_equal(mpn(positive = c(5, 1), tubes = tubes, amount = amount)$MPN,
               3.48, tol = .01)
})

test_that("Blodgett - RI matches Table 4, Blodgett (2002)", {
  tubes  <- c(5, 5, 5)
  amount <- c(1, .1, .01)
  expect_equal(mpn(positive = c(3, 0, 0), tubes = tubes, amount = amount)$RI,
               1, tol = .001)
  expect_equal(mpn(positive = c(2, 1, 0), tubes = tubes, amount = amount)$RI,
               .354, tol = .01)
  expect_equal(mpn(positive = c(0, 4, 1), tubes = tubes, amount = amount)$RI,
               5.85e-07, tol = 1e-07)
  expect_equal(mpn(positive = c(3, 4, 0), tubes = tubes, amount = amount)$RI,
               2.51e-03, tol = 1e-04)
  expect_equal(mpn(positive = c(3, 3, 1), tubes = tubes, amount = amount)$RI,
               2.34e-03, tol = 1e-04)
  expect_equal(mpn(positive = c(0, 2, 5), tubes = tubes, amount = amount)$RI,
               2.91e-13, tol = 1e-14)
})


# Jarvis (2010) ----------------------------------------------------------------

test_that("Jarvis - some positive, some negative", {
  #Table 1
  x <- mpn(positive = c(2, 2, 0), tubes = c(3, 3, 3), amount = c(1, .1, .01))
  expect_equal(x$MPN, 2.1, tol = .1)
  expect_equal(x$LB, .71, tol = .025)
  expect_equal(x$UB, 6.2, tol = .1)
  expect_equal(x$RI, .069, tol = .001)

  y <- mpn(positive = c(3, 1, 1), tubes = c(3, 3, 3), amount = c(1, .1, .01))
  expect_equal(y$MPN, 7.5, tol = .1)
  expect_equal(y$LB, 1.9, tol = .1)
  expect_equal(y$UB, 30, tol = 1)
  expect_equal(y$RI, .209, tol = .001)

  #Table 2
  z <- mpn(positive = c(50, 34, 7), tubes = c(50, 49, 49),
           amount = 2 * c(1, .1, .01))
  expect_equal(z$MPN, 6.2, tol = .1)
  expect_equal(z$LB, 4.5, tol = .1)
  expect_equal(z$UB, 8.6, tol = .1)
  expect_equal(z$RI, .746, tol = .001)
})

test_that("Jarvis - all negative tubes", {
  #Jarvis Table 1
  x <- mpn(positive = c(0, 0, 0), tubes = c(3, 3, 3), amount = c(1, .1 , .01),
           conf_level = .975)  #Jarvis used alpha / 2
  expect_equal(x$MPN, 0)
  expect_equal(x$variance, NA)
  expect_equal(x$var_log, NA)
  expect_equal(x$LB, 0)
  expect_equal(x$UB, 1.1, tol = .1)
  expect_equal(x$RI, 1)
})

test_that("Jarvis - all positive tubes", {
  #Table 1
  x <- mpn(positive = c(3, 3, 3), tubes = c(3, 3, 3), amount = c(1, .1, .01),
           conf_level = .975)  #Jarvis used alpha / 2
  expect_true(is.infinite(x$MPN))
  expect_true(is.na(x$variance))
  expect_true(is.na(x$var_log))
  expect_equal(x$LB, 36, tol = .05)
  expect_true(is.infinite(x$UB))
  expect_equal(x$RI, 1)
})


# BAM tables -------------------------------------------------------------------

test_that("BAM tables - some positive, some negative", {
  amount <- c(.1, .01, .001)
  #Table 1
  tubes1 <- c(3, 3, 3)
  x1 <- mpn(positive = c(0, 0, 1), tubes = tubes1, amount = amount)
  expect_equal(signif(x1$MPN, 2), 3.0, tol = .001)
  y1 <- mpn(positive = c(1, 1, 1), tubes = tubes1, amount = amount)
  expect_equal(signif(y1$MPN, 2), 11, tol = .001)
  z1 <- mpn(positive = c(2, 2, 0), tubes = tubes1, amount = amount)
  expect_equal(signif(z1$MPN, 2), 21, tol = .001)
  q1 <- mpn(positive = c(3, 3, 2), tubes = tubes1, amount = amount)
  expect_equal(signif(q1$MPN, 2), 1100, tol = .001)

  #Table 2
  tubes2 <- c(5, 5, 5)
  x2 <- mpn(positive = c(0, 0, 1), tubes = tubes2, amount = amount)
  expect_equal(signif(x2$MPN, 2), 1.8, tol = .001)
  y2 <- mpn(positive = c(1, 2, 0), tubes = tubes2, amount = amount)
  expect_equal(signif(y2$MPN, 2), 6.1, tol = .001)
  z2 <- mpn(positive = c(4, 0, 2), tubes = tubes2, amount = amount)
  expect_equal(signif(z2$MPN, 2), 21, tol = .001)
  q2 <- mpn(positive = c(5, 5, 4), tubes = tubes2, amount = amount)
  expect_equal(signif(q2$MPN, 2), 1600, tol = .001)

  #Table 3
  tubes3 <- c(10, 10, 10)
  x3 <- mpn(positive = c(0, 0, 1), tubes = tubes3, amount = amount)
  expect_equal(signif(x3$MPN, 2), 0.9, tol = .001)
  y3 <- mpn(positive = c(3, 2, 2), tubes = tubes3, amount = amount)
  expect_equal(signif(y3$MPN, 2), 7.5, tol = .001)
  z3 <- mpn(positive = c(8, 2, 0), tubes = tubes3, amount = amount)
  expect_equal(signif(z3$MPN, 2), 17, tol = .001)
  q3 <- mpn(positive = c(10, 10, 9), tubes = tubes3, amount = amount)
  expect_equal(signif(q3$MPN, 2), 2300, tol = .001)

  #Table 5
  x4 <- mpn(positive = 1, tubes = 10, amount = 10)
  expect_equal(signif(100 * x4$MPN, 2), 1.1, tol = .01)
  y4 <- mpn(positive = 5, tubes = 10, amount = 10)
  expect_equal(signif(100 * y4$MPN, 2), 6.9, tol = .01)
  z4 <- mpn(positive = 9, tubes = 10, amount = 10)
  expect_equal(signif(100 * z4$MPN, 2), 23, tol = .1)
})
