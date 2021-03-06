#' Calculate aerobic plate count (APC)
#'
#' \code{apc} calculates the Aerobic Plate Count (\emph{APC}) point estimate
#' and confidence interval of colony forming units (CFU). Adjusts for
#' too-numerous-to-count (TNTC) plates using the maximum likelihood method of
#' Haas et al. (2014).
#'
#' @param count A vector of CFU counts in each scorable (countable) plate.
#' @param amount_scor A vector of inoculum amounts (in ml) in each scorable
#'   plate. See \emph{Details} section.
#' @param amount_tntc A vector of inoculum amounts (in ml) in each TNTC plate.
#' @param tntc_limit A vector (or scalar) of the limit above which the plate
#'   counts are considered too-numerous-to-count (often 100, 250, or 300). Each
#'   plate can potentially have a different value. Default is 250.
#' @param conf_level A scalar value between zero and one for the confidence
#'   level. Typically 0.95 (i.e., a 95 percent confidence interval).
#'
#' @return A list containing:
#'   \itemize{
#'     \item{\strong{APC}: }{The aerobic plate count point estimate in CFU/ml.}
#'     \item{\strong{conf_level}: }{The confidence level used.}
#'     \item{\strong{LB}: }{The lower bound of the confidence interval.}
#'     \item{\strong{UB}: }{The upper bound of the confidence interval.}
#'   }
#'
#' @details As an example, assume we start with four plates and 1 ml of
#'   undiluted inoculum. For the first two plates we use a 100-fold dilution;
#'   for the other two plates we use a 1,000-fold dilution. The first two plates
#'   were TNTC with limits of 300 and 250. The other plates had CFU counts of
#'   28 and 20. We now have
#'   \code{count = c(28, 20)},
#'   \code{amount_scor = 1 * c(.001, .001)},
#'   \code{amount_tntc = 1 * c(.01, .01)}, and
#'   \code{tntc_limit = c(300, 250)}.
#'
#' @details Currently, confidence intervals can only be calculated using the
#'   likelihood ratio (LR) approach described in Haas et al. (2014).
#'
#' @section Warnings:
#'   The likelihood ratio confidence interval assumptions depend on asymptotic
#'   theory. Therefore, the confidence interval results will be better with
#'   larger experiments.
#'
#' @section Warnings:
#'   \code{apc()} will fail in certain cases where the TNTC results are
#'   extremely unlikely to occur when taking the scorable (countable) plates
#'   into consideration. In other words, if the countable plates suggest a low
#'   concentration of microbes, then TNTC plates at higher dilution levels are
#'   probably due to experimental error. Mathematically, the probability is so
#'   small that the likelihood function is essentially zero.
#'
#' @examples
#' #------- "Quantitative Microbial Risk Assessment (Haas et al., 2014) --------
#'
#' # Table 6.1 (Sample A)
#' my_count <- c(1, 2, 1, 0, 0, 1, 1, 3, 6, 8, 4)
#' my_amount_scor <- c(1, 1, 1, 1, 1, 2.5, 2.5, 2.5, 2.5, 5, 5)
#' apc(my_count, my_amount_scor)  #1.08
#'
#' # Table 6.1 (Sample B)
#' my_count <- c(1, 0, 5, 1, 0, 5, 0, 1, 5, 1, 8)
#' my_amount_scor <- c(1, 1, 1, 1, 1, 2.5, 2.5, 2.5, 2.5, 5, 5)
#' apc(my_count, my_amount_scor)  #1.08
#'
#' # Table 6.2
#' my_count <- c(12, 8, 15, 40, 58)
#' my_amount_scor <- c(1, 1, 1, 10, 10)
#' my_amount_tntc <- c(10, 100, 100, 100)
#' my_tntc_limit <- 100
#' apc(my_count, my_amount_scor, my_amount_tntc, my_tntc_limit) #~7 (6.03, 7.96)
#'
#'
#' #----------- "Averaging of TNTC Counts" (Haas & Heller, 1988) ---------------
#' # Note:
#' #  Results are slightly different due mostly to differences in how the TNTC
#' #  portion of the likelihood function is formulated (i.e., incomplete gamma
#' #  function vs. infinite Poisson sum--see Haas et al. (2014) for details of
#' #  this mathematical relationship).
#'
#' my_count <- c(10, 12, 23, 48, 63)
#' my_amount_scor <- c(1, 1, 1, 5, 5)
#' my_amount_tntc <- c(5, 10, 10)
#' my_tntc_limit <- 80
#' apc(my_count, my_amount_scor, my_amount_tntc, my_tntc_limit)
#' #Haas & Heller: APC = 13.28 CFU/ml
#'
#' @references Bacteriological Analytical Manual, 8th Edition, Chapter 3,
#'   \url{https://www.fda.gov/food/foodscienceresearch/laboratorymethods/ucm063346.htm}
#'
#' @references Haas CN, Heller B (1988). "Averaging of TNTC counts."
#'   \emph{Applied and Environmental Microbiology}, 54(8), 2069-2072.
#'   \url{https://aem.asm.org/content/54/8/2069}
#'
#' @references Haas CN, Rose JB, Gerba CP (2014). "Quantitative microbial risk
#'   assessment, Second Ed." \emph{John Wiley & Sons, Inc.},
#'   ISBN 978-1-118-14529-6.
#'
#' @seealso \code{\link{mpn}} for Most Probable Number
#' @importFrom stats uniroot
#' @importFrom stats optimize
#' @importFrom stats pgamma
#' @importFrom stats qchisq
#' @export

apc <- function(count, amount_scor, amount_tntc = NULL, tntc_limit = 250,
                conf_level = 0.95) {

  .checkInputs_apc(count, amount_scor, amount_tntc, tntc_limit, conf_level)

  APC <- .ptEst_APC(count, amount_scor, amount_tntc, tntc_limit)

  like_ratio_CI <- .likeRatioCI_APC(lambda_hat = APC,
                                    count = count,
                                    amount_scor = amount_scor,
                                    amount_tntc = amount_tntc,
                                    tntc_limit = tntc_limit,
                                    conf_level = conf_level)
  LB <- like_ratio_CI$LB
  UB <- like_ratio_CI$UB

  list(APC = APC, conf_level = conf_level, LB = LB, UB = UB)
}
