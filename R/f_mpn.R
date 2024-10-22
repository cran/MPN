#' Calculate most probable number (MPN)
#'
#' \code{mpn} calculates the Most Probable Number (\emph{MPN}) point estimate
#' and confidence interval for microbial concentrations. Also calculates
#' Blodgett's (2002, 2005, 2010) Rarity Index (\emph{RI}).
#'
#' @param positive A vector of number of positive tubes at each dilution level.
#' @param tubes A vector of total number of tubes at each dilution level.
#' @param amount A vector of the amount of inoculum per tube at each dilution
#'   level. See \emph{Details} section.
#' @param conf_level A scalar value between zero and one for the confidence
#'   level. Typically 0.95 (i.e., a 95 percent confidence interval).
#' @param CI_method The method used for calculating the confidence interval.
#'   Choices are \code{"Jarvis"} or \code{"LR"} (likelihood ratio). See
#'   \emph{Details} section.
#' @param tol A scalar value for tolerance to be passed to
#'   \code{stats::uniroot()}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \strong{MPN}: The most probable number point estimate for
#'       microbial density (concentration).
#'     \item \strong{MPN_adj}: The bias-adjusted point estimate for MPN.
#'     \item \strong{variance}: The estimated variance (see Jarvis et al.) of
#'       the \emph{MPN} estimate if \code{CI_method = "Jarvis"}. If all tubes
#'       are positive or all negative, \code{variance} will be \code{NA}. If
#'       \code{CI_method} is not \code{"Jarvis"}, \code{variance} will be
#'       \code{NA}.
#'     \item \strong{var_log}: The estimated variance of the natural log of
#'       the MPN estimate (see Jarvis et al.) using the Delta Method. If all
#'       tubes are positive or all negative, \code{var_log} will be \code{NA}.
#'       If \code{CI_method} is not \code{"Jarvis"}, \code{var_log} will be
#'       \code{NA}.
#'     \item \strong{conf_level}: The confidence level used.
#'     \item \strong{CI_method}: The confidence interval method used.
#'     \item \strong{LB}: The lower bound of the confidence interval.
#'     \item \strong{UB}: The upper bound of the confidence interval.
#'     \item \strong{RI}: The rarity index.
#'   }
#'
#' @details As an example, assume we start with 3g of undiluted inoculum per
#'   tube, then use a 10-fold dilution for 2 dilutions. We now have
#'   \code{amount = 3 * c(1, .1, .01)}.
#'
#' @details When all tubes are negative, the point estimate of \emph{MPN} is
#'   zero (same approach as Jarvis et al.). The point estimate for the BAM
#'   tables "is listed as less than the lowest MPN for an outcome with at least
#'   one positive tube" (App.2).
#'
#' @details When all tubes are positive, the point estimate for \emph{MPN} is
#'   \code{Inf} (same approach as Jarvis et al.) since no finite maximum
#'   likelihood estimate (MLE) exists. The BAM tables "list the MPN for this
#'   outcome as greater than the highest MPN for an outcome with at least one
#'   negative tube" (App.2). Here, the difference is probably trivial since the
#'   sample should be further diluted if all tubes test positive.
#'
#' @details The bias adjustment for the point estimate uses the method of Salama
#'   et al. (1978). Also see Haas (1989).
#'
#' @details Currently, confidence intervals can only be calculated using the
#'   Jarvis (2010) or likelihood ratio (LR) approach (Ridout, 1994). The BAM
#'   tables use an alternate approach. We slightly modified Jarvis' approach
#'   when all tubes are positive or all are negative; we use \eqn{\alpha}
#'   instead of \eqn{\alpha / 2} since these are one-sided intervals. The Ridout
#'   (1994) LR approach uses the same technique (with \eqn{\alpha}) for these
#'   two extreme cases.
#'
#' @details If the Rarity Index is less than \code{1e-04}, the experimental
#'   results are highly improbable. The researcher may consider running the
#'   experiment again and/or changing the dilution levels.
#'
#' @section Warnings:
#'   The Jarvis confidence interval assumptions of approximate normality (Delta
#'   Method and asymptotic normality of maximum likelihood estimators) depend on
#'   large-sample theory. The likelihood ratio assumptions also depend on
#'   large-sample theory. Therefore, the Jarvis and LR confidence interval
#'   approaches work best with larger experiments.
#'
#' @examples
#' # Compare MPN, 95% CI, and RI to Jarvis -------------------------------------
#'
#' # Table 1
#' mpn(positive = c(3, 1, 1), tubes = c(3, 3, 3), amount = c(1, .1, .01))
#'   #Jarvis: 7.5 (1.9, 30) RI = .209
#'
#' mpn(positive = c(0, 0, 0), tubes = c(3, 3, 3), amount = c(1, .1, .01))
#'   #Jarvis: 0 (0, 1.1) RI = 1
#' mpn(positive = c(0, 0, 0), tubes = c(3, 3, 3), amount = c(1, .1, .01),
#'     conf_level = .975)$UB  #alpha / 2
#'
#' mpn(positive = c(3, 3, 3), tubes = c(3, 3, 3), amount = c(1, .1, .01))
#'   #Jarvis: Inf (36, Inf) RI = 1
#' mpn(positive = c(3, 3, 3), tubes = c(3, 3, 3), amount = c(1, .1, .01),
#'     conf_level = .975)$LB  #alpha / 2
#'
#' # Table 2
#' mpn(positive = c(20, 14, 3), tubes = c(20, 20, 20), amount = c(1, .1, .01))
#'   #Jarvis: 13 (7.6, 21) RI = 0.794
#'
#' mpn(positive = c(50, 35, 7), tubes = c(50, 50, 50),
#'     amount = 2 * c(1, .1, .01))
#'   #Jarvis: 6.3 (4.5, 8.7) RI = .806
#'
#' mpn(positive = c(1, 5, 3, 1, 1), tubes = c(1, 5, 5, 5, 5),
#'     amount = c(5, 1, .5, .1, .05))
#'   #Jarvis: 2.7 (1.3, 5.5) RI = .512
#'
#' # Compare MPN and 95% CI to BAM tables --------------------------------------
#'
#' # Table 1
#' mpn(positive = c(0, 0, 0), tubes = c(3, 3, 3), amount = c(.1, .01, .001))
#'   #BAM: <3.0 (-, 9.5)
#'
#' mpn(positive = c(0, 0, 1), tubes = c(3, 3, 3), amount = c(.1, .01, .001))
#'   #BAM: 3.0 (0.15, 9.6)
#'
#' mpn(positive = c(2, 2, 0), tubes = c(3, 3, 3), amount = c(.1, .01, .001))
#'   #BAM: 21 (4.5, 42)
#'
#' mpn(positive = c(3, 3, 3), tubes = c(3, 3, 3), amount = c(.1, .01, .001))
#'   #BAM: >1100 (420, -)
#' mpn(positive = c(3, 3, 2), tubes = c(3, 3, 3), amount = c(.1, .01, .001))$MPN
#'
#'
#' # Table 2
#' mpn(positive = c(0, 0, 0), tubes = c(5, 5, 5), amount = c(.1, .01, .001))
#'   #BAM: <1.8 (-, 6.8)
#' mpn(positive = c(0, 0, 1), tubes = c(5, 5, 5), amount = c(.1, .01, .001))$MPN
#'
#' mpn(positive = c(4, 0, 2), tubes = c(5, 5, 5), amount = c(.1, .01, .001))
#'   #BAM: 21 (6.8, 40)
#'
#' mpn(positive = c(5, 5, 5), tubes = c(5, 5, 5), amount = c(.1, .01, .001))
#'   #BAM: >1600 (700, -)
#' mpn(positive = c(5, 5, 4), tubes = c(5, 5, 5), amount = c(.1, .01, .001))$MPN
#'
#' # Compare MPN and 95% LR CI to Ridout (1994) --------------------------------
#'
#' # Table 1
#' mpn(positive = c(0, 0, 0), tubes = c(3, 3, 3), amount = c(.1, .01, .001),
#'     CI_method = "LR")
#'   #Ridout: 0 (0, 9.0)
#' mpn(positive = c(2, 2, 0), tubes = c(3, 3, 3), amount = c(.1, .01, .001),
#'     CI_method = "LR")
#'   #Ridout: 21.1 (6.2, 54.3)
#' mpn(positive = c(3, 3, 3), tubes = c(3, 3, 3), amount = c(.1, .01, .001),
#'     CI_method = "LR")
#'   #Ridout: Inf (465.1, Inf)
#'
#' @references Bacteriological Analytical Manual, 8th Edition, Appendix 2,
#'   \url{https://www.fda.gov/food/laboratory-methods-food/bam-appendix-2-most-probable-number-serial-dilutions}
#'
#' @references Blodgett RJ (2002). "Measuring improbability of outcomes from a
#'   serial dilution test." \emph{Communications in Statistics: Theory and
#'   Methods}, 31(12), 2209-2223.
#'
#' @references Blodgett RJ (2005). "Serial dilution with a confirmation step."
#'   \emph{Food Microbiology}, 22(6), 547-552.
#'
#' @references Blodgett RJ (2010). "Does a serial dilution experiment's model
#'   agree with its outcome?" \emph{Model Assisted Statistics and Applications},
#'   5(3), 209-215.
#'
#' @references Haas CN (1989). "Estimation of microbial densities from dilution
#'   count experiments" \emph{Applied and Environmental Microbiology} 55(8),
#'   1934-1942.
#'
#' @references Haas CN, Rose JB, Gerba CP (2014). "Quantitative microbial risk
#'   assessment, Second Ed." \emph{John Wiley & Sons, Inc.},
#'   ISBN 978-1-118-14529-6.
#'
#' @references Jarvis B, Wilrich C, Wilrich P-T (2010). "Reconsideration of the
#'   derivation of Most Probable Numbers, their standard deviations, confidence
#'   bounds and rarity values." \emph{Journal of Applied Microbiology}, 109,
#'   1660-1667.
#'
#' @references Ridout MS (1994). "A comparison of confidence interval methods
#'   for dilution series experiments." \emph{Biometrics}, 50(1), 289-296.
#'
#' @references Salama IA, Koch GG, Tolley DH. (1978) "On the estimation of the
#'   most probable number in a serial dilution technique." \emph{Communications
#'   in Statistics - Theory and Methods}, 7(13), 1267-1281.
#'
#' @seealso Shiny app:  \url{https://pub-connect.foodsafetyrisk.org/microbial/mpncalc/}
#' @seealso \code{\link{apc}} for Aerobic Plate Count
#' @importFrom stats uniroot
#' @importFrom stats dbinom
#' @importFrom stats qnorm
#' @importFrom stats qchisq
#' @export

mpn <- function(positive, tubes, amount, conf_level = 0.95,
                CI_method = c("Jarvis", "LR"), tol = 1e-06) {

  .checkInputs_mpn(positive = positive, tubes = tubes, amount = amount,
                   conf_level = conf_level, tol = tol)
  CI_method <- match.arg(CI_method)

  MPN <- .ptEst_MPN(positive, tubes, amount, tol)
  MPN_adj <- .ptEstAdj_MPN(MPN, tubes, amount)

  variance   <- NA
  var_logMPN <- NA
  if (CI_method == "Jarvis") {
    jarvis     <- .jarvisCI_MPN(MPN = MPN, positive = positive, tubes = tubes,
                                amount = amount, conf_level = conf_level,
                                tol = tol)
    variance   <- jarvis$variance
    var_logMPN <- jarvis$var_logMPN
    LB <- jarvis$LB
    UB <- jarvis$UB
  } else {
    like_ratio <- .likeRatioCI_MPN(MPN = MPN, positive = positive,
                                   tubes = tubes, amount = amount,
                                   conf_level = conf_level, tol = tol)
    LB <- like_ratio$LB
    UB <- like_ratio$UB
  }

  rarity <- .rarity_MPN(MPN, positive, tubes, amount)

  list(MPN = MPN, MPN_adj = MPN_adj, variance = variance, var_log = var_logMPN,
       conf_level = conf_level, CI_method = CI_method, LB = LB, UB = UB,
       RI = rarity)
}
