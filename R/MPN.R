#' \pkg{MPN}: Most Probable Number and Other Microbial Enumeration Techniques
#'
#' \pkg{MPN} is a package for calculating the Most Probable Number (MPN) to
#' quantify the concentration of microbes in serial dilutions of laboratory
#' samples. Also calculates the Aerobic Plate Count (APC) for similar
#' experiments.
#'
#' @section Functions:
#' The function \code{mpn} calculates the Most Probable Number (\emph{MPN})
#' point estimate and confidence interval for microbial concentrations. Also
#' calculates Blodgett's (2002, 2005, 2010) Rarity Index (\emph{RI}). The MPN
#' calculation is described in the Bacteriological Analytical Manual (BAM, 8th
#' ed., Appendix 2) and Jarvis et al. (2010).
#'
#' @section Functions:
#' \code{apc} calculates the Aerobic Plate Count (\emph{APC}) point estimate
#' and confidence interval of colony forming units (CFU). Adjusts for
#' too-numerous-to-count (TNTC) plates using the maximum likelihood method of
#' Haas et al. (2014).
#'
#' @references Bacteriological Analytical Manual, 8th Edition, Appendix 2,
#'   \url{https://www.fda.gov/food/laboratory-methods-food/bam-appendix-2-most-probable-number-serial-dilutions}
#'
#' @references Bacteriological Analytical Manual, 8th Edition, Chapter 3,
#'   \url{https://www.fda.gov/food/laboratory-methods-food/bam-chapter-3-aerobic-plate-count}
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
#' @references Haas CN, Heller B (1988). "Averaging of TNTC counts."
#'   \emph{Applied and Environmental Microbiology}, 54(8), 2069-2072.
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
#' @keywords internal
"_PACKAGE"

NULL
