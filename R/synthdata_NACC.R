#' synthdata NACC
#'
#' @importFrom psych ICC
#' @name synthdata_NACC
#' @docType data
#'
#' @description NACC MoCA synthetic example data (2433 observations of patients
#'   with no clinical assessment of cognitive impairment and 2644 observations
#'   with a clinical assessment of some form of cognitive impairment.
#'
#'   \itemize{ \item ID patient id; sequential, randomly assigned \item center
#'   Alphanumeric id of the clinical center where the data has been collected
#'   (30 centers) \item ref.1 Gold standard (true status) at the first
#'   measurement. 0: no cognitive impairment 1: cognitive impairment \item
#'   MOCATOTS.1 Total MoCA score at the first measurement (0 .. 30) \item
#'   vdate.1 Date of the first measurement \item ref.2 Gold standard (true
#'   status) at the second measurement. 0: no cognitive impairment 1: cognitive
#'   impairment \item MOCATOTS.2 Total MoCA score at the second measurement (0
#'   .. 30) \item vdate.2 Date of the second measurement }
#'
#' @details For use as an example, a single data set of 6670 observations is
#'   generated based on the NACC dataset, from 30 different clinical centers. To
#'   generate the artificial data, the R package synthpop (Nowok B, Raab GM,
#'   Dibben C, 2016) is used to create the synthetic data, based on the original
#'   data from the Uniform Data Set (UDS), collected by the University of
#'   Washington’s National Alzheimer’s Coordinating Center (NACC). The syntetic
#'   data provide similar statistical results, but differ for each individual
#'   and each clinical center.
#'   These data are provided as data for the replication of the examples.
#'   Results of the real data are presented in Landsheer (In Press).
#'
#'   Researchers who want to use these data for other purposes than replication
#'   of the results presented here, are kindly requested to submit a new request
#'   for the original data to the NACC. The user of the data may either get a
#'   new file or request a file using the specifications of the original data
#'   file (https://www.alz.washington.edu/).
#'
#' @references Nowok B, Raab GM, Dibben C (2016). “synthpop: Bespoke Creation of
#'   Synthetic Data in R.” Journal of Statistical Software, 74(11), 1–26.
#'   doi:10.18637/jss.v074.i11.
#'
#'   Landsheer, J. A. (In press). Impact of the Prevalence of Cognitive
#'   Impairment on the Accuracy of the Montreal Cognitive Assessment: The
#'   advantage of using two MoCA thresholds to identify error-prone test scores.
#'   Alzheimer Disease and Associated Disorders.
#'   https://doi.org/10.1097/WAD.0000000000000365
#'
#' @examples
#' data(synthdata_NACC) # needs R version 3.5 or later
#' head(synthdata_NACC) # Show head of the dataset
#' nrow(synthdata_NACC) # total number of observations
#' # select part of data for the first measurement
#' # N.B. ref is not available when it is inconclusive
#' m1 = synthdata_NACC[!is.na(synthdata_NACC$MOCATOTS.1)
#'                    & !is.na(synthdata_NACC$ref.1), ]
#' # preliminary check data for possible missing values
#' addmargins(table(m1$ref.1, m1$MOCATOTS.1, useNA = 'always'))
#' # Show the data
#' barplotMD(m1$ref.1, m1$MOCATOTS.1)
#'
#' # calculate the difference between the two measurements in days
#' ddiff = (m1$vdate.2 - m1$vdate.1)
#' # There is a wide variety !!!
#' summary(ddiff)
#' # Estimate the test-retest reliability
#' library(psych)
#' ICC(na.omit(cbind(m1$MOCATOTS.1, m1$MOCATOTS.2)))
#' # Reducing the variety of time between measurements:
#' timesel = (ddiff >= 335) & (ddiff <= 395)
#' ICC(na.omit(cbind(m1$MOCATOTS.1[timesel], m1$MOCATOTS.2[timesel])))
#'
#' # error when using default calculated value for roll.length
#' # RPV(m1$ref.1, m1$MOCATOTS.1, reliability = .86)
#' RPV(m1$ref.1, m1$MOCATOTS.1, reliability = .86, roll.length = 5)

#'
#' @keywords data
NULL
