#' Primitive non-parametric function for calculating AUC or C-statistic for two
#' comparable samples with ordinal data.
#'
#' @param norm Ordinal data of the norm group (controls).
#' @param abnorm Ordinal data of the abnorm group (patients).
#'
#' @return The statistic AUC (Area under the Receiver Operating Characteristics
#'   Curve), also known as the C-statistic or Concordance statistic.
#' @details This function does not check anything. Argument \code{norm} and
#'   \code{abnorm} must have the correct ordinal data. When \code{mean(norm) >
#'   mean(abnorm)} it is assumed that lower sores indicate deviation from the
#'   norm, and otherwise that higher scores indicate deviation of the norm.
#'   This function can handle very large files.
#'
#' @export
#'
#' @examples
#' norm = round(rnorm(100, 3, 1))
#' abnorm= round(rnorm(80, 5, 2))
#' simple_auc(norm, abnorm)

simple_auc <- function(norm, abnorm){
  if(length(norm)==0 | length(abnorm)==0) return(NA)
  n1 = length(abnorm);
  n0 = length(norm);
  r = rank(c(abnorm, norm))
  auc = (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n0)
  if (mean(norm) > mean(abnorm))
    return(1-auc) 
  else return(auc)
}
