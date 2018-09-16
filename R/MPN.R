#' \pkg{MPN}: Most Probable Number for Serial Dilutions
#'
#' \pkg{MPN} is a package for calculating the Most Probable Number (MPN) for
#' serial dilutions to quantify the concentration of microbes in a laboratory
#' sample.
#'
#' @section Functions:
#' The function \code{\link{mpn}} calculates the Most Probable Number
#' (\emph{MPN}). The MPN calculation is described in the Bacteriological
#' Analytical Manual (BAM, 8th ed., Appendix 2) and Jarvis et al. (2010).
#' Calculates variance, variance of the natural log, and confidence interval
#' using the approaches described in Jarvis et al. (2010). Also calculates
#' Blodgett's (2002, 2005, 2010) Rarity Index (\emph{RI}).
#'
#' @references Bacteriological Analytical Manual, 8th Edition, Appendix 2,
#'   \url{https://www.fda.gov/Food/FoodScienceResearch/LaboratoryMethods/ucm109656.htm}
#'
#' @references Blodgett RJ (2002). "Measuring improbability of outcomes from a
#'   serial dilution test." \emph{Communications in Statistics: Theory and
#'   Methods}, 31(12), 2209-2223. \url{https://doi.org/10.1081/STA-120017222}
#'
#' @references Blodgett RJ (2005). "Serial dilution with a confirmation step."
#'   \emph{Food Microbiology}, 22(6), 547-552.
#'   \url{https://doi.org/10.1016/j.fm.2004.11.017}
#'
#' @references Blodgett RJ (2010). "Does a serial dilution experiment's model
#'   agree with its outcome?" \emph{Model Assisted Statistics and Applications},
#'   5(3), 209-215. \url{https://doi.org/10.3233/MAS-2010-0157}
#'
#' @references Jarvis B, Wilrich C, Wilrich P-T (2010). "Reconsideration of the
#'   derivation of Most Probable Numbers, their standard deviations, confidence
#'   bounds and rarity values." \emph{Journal of Applied Microbiology}, 109,
#'   1660-1667. \url{https://doi.org/10.1111/j.1365-2672.2010.04792.x}
#'
#' @docType package
#' @name MPN
NULL