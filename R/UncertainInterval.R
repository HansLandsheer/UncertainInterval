#' @title Set of functions for the determination of an Uncertain Interval of
#'   test scores
#'
#' @description A collection of functions to determine a range of test scores
#'   that are inconclusive and do not allow a diagnosis (other than Uncertain)
#'   and to access its qualities.
#'
#' @seealso \code{\link{ui.nonpar}},  \code{\link{plotMD}},
#'   \code{\link{get.intersection}},  \code{\link{quality.threshold}},
#'   \code{\link{quality.threshold.uncertain}}
#' @name UncertainInterval
#' @details Uncertain test scores are scores that have about the same density in
#'   the two distributions of patients with and without the targeted condition.
#'   This range is typically found around the optimal cut-point, that is the
#'   point of intersection or Youden index (Schisterman et al., 2005).
#'
#'   \describe{ \item{Criteria}{Most functions in this package use a specified
#'   low value for the sensitivity and specificity of the test scores within the
#'   uncertain interval to find this uncertain interval (default UI.Se = UI.Sp =
#'   .55). The most recent added function \code{\link{RPV}} for ordinal test
#'   scores uses the odds of the target condition of near 1 to identify the
#'   uncertain interval (default < 2). This library also contains two
#'   alternative definitions. 1. Coste et al. (2003) defined a grey zone in
#'   between positive and negative conclusions (see \code{\link{greyzone}}),
#'   minimum desired values for respectively the positive and negative post-test
#'   probability, with defaults .95 and .05. 2. Greiner (1995) defined a middle
#'   inconclusive zone of intermediate values (see \code{\link{TG.ROC}}), with
#'   desired minimum values for dichotomous Se and Sp, with default values of
#'   .9. See Index for all available functions and plot possibilities.}
#'   \item{Glossary}{ In general, the prefix MCI is used when a statistic is
#'   calculated for the test scores that are used for a positive or negative
#'   classification. The prefix UI is used when the statistic is applied to the
#'   test scores in the uncertain interval.   
#'   \itemize{ \item{Se and Sp}{Se and Sp are statistics
#'   that are developed for a single dichotomous cut-point.} 
#'   \item{MCI.Se and
#'   MCI.Sp}{ Sensitivity and specificity calculated for the More
#'   Certain Intervals (MCIs) outside the Uncertain Interval (UI), that is,
#'   omitting the test scores in the UI. The meaning of Se and Sp changes from sensitivity
#'   and specificity of the test (or all test scores) to sensitivy and specificity of the test scores
#'   used for classification.} 
#'   \item{UI.Se and
#'   UI.Sp}{Sensitivity and specificity for the test scores inside the uncertain
#'   interval. Please note that the uncertain interval always falls around the
#'   point of intersection (optimal threshold or Youden threshold) and that for
#'   the calculation of UI.Se and UI.Sp the point of intersection is used as
#'   threshold within the uncertain interval.} 
#'   \item{NPV and PPV}{Predictive values for respectively the negative and the
#'   positive class. Can be used with both dichotomous and trichotomous sections
#'   of the test scores. The prefix MCI is sometimes used, but is superfluous.}
#'   \item{PV.class}{Predictive value for
#'   class when the meaning of class is selfexplanatory.}
#'   \item{NPV.class and PPV.class} {Negative and Positive Predictive value when the scores in class
#'   are used for a negative, respectively positive classification. When
#'   predictive values are calculated for the same class, NPV.class = 1 -
#'   PPV.class. }
#'   \item{NPV.ui and PPV.ui} {Negative and Positive Predictive value when all test
#'   scores in the uncertain interval would be used for a negative, respectively
#'   positive classification. These values can be expected to be close to .5.
#'   When predictive values are calculated for the same class, NPV.ui = 1 -
#'   PPV.ui.} 
#'   \item{UI.NPV and UI.PPV}{Negative and Positive Predictive value when test
#'   scores in the uncertain interval respectively above and below the point of
#'   intersection would be used for a negative, respectively positive
#'   classification. These values can be expected to be close to .5, but 
#'   slightly higher than NPV.ui and PPV.ui}
#'   \item{SNPV, SPPV, SPV.class, SNPV.class, SPPV.class, SNPV.ui, SPPV.ui,
#'   UI.SNPV and UI.SPPV}{The standardized versions of the predictive values
#'   mentioned above.}
#'   } } }
#' @references Landsheer, J. A. (2016). Interval of Uncertainty: An Alternative
#'   Approach for the Determination of Decision Thresholds, with an Illustrative
#'   Application for the Prediction of Prostate Cancer. PloS One, 11(11),
#'   e0166007.
#'
#'   Landsheer, J. A. (2018). The Clinical Relevance of Methods for Handling
#'   Inconclusive Medical Test Results: Quantification of Uncertainty in Medical
#'   Decision-Making and Screening. Diagnostics, 8(2), 32.
#'   https://doi.org/10.3390/diagnostics8020032
#'
#'   Schisterman, E. F., Perkins, N. J., Liu, A., & Bondell, H. (2005). Optimal
#'   cut-point and its corresponding Youden Index to discriminate individuals
#'   using pooled blood samples. Epidemiology, 73-81.
#'
#'   Greiner, M. (1995). Two-graph receiver operating characteristic (TG-ROC): A
#'   Microsoft-EXCEL template for the selection of cut-off values in diagnostic
#'   tests. Journal of Immunological Methods, 185(1), 145-146.
#'
#'   Coste, J., & Pouchot, J. (2003). A grey zone for quantitative diagnostic
#'   and screening tests. International Journal of Epidemiology, 32(2), 304-313.


NULL
#> NULL
## "_PACKAGE"
##> [1] "_PACKAGE"
