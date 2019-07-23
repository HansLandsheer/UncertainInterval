#'Function for the determination of an inconclusive interval for continuous test
#'scores
#'
#'@name ui.nonpar
#'@description This function uses a non-parametric approach to determine an
#'  interval around the intersection of the two distributions of individuals
#'  without (0) and with (1) the targeted condition. The Uncertain Interval is
#'  generally defined as an interval below and above the intersection, where the
#'  densities of the two distributions of patients with and without the targeted
#'  condition are about equal. These test scores are considered as inconclusive
#'  for the decision for or against the targeted condition. The interval is
#'  restricted both by a maximum sensitivity of the test scores within the
#'  uncertain interval (sens.ui) and by a maximum specificity of the test scores
#'  within the uncertain interval (spec.ui).
#'@param ref The reference standard. A column in a data frame or a vector
#'  indicating the classification by the reference test. The reference standard
#'  must be coded either as 0 (absence of the condition) or 1 (presence of the
#'  condition)
#'@param test The index test or test under evaluation. A column in a dataset or
#'  vector indicating the test results in a continuous scale.
#'@param sens.ui (default = .55). The sensitivity of the test scores within the
#'  uncertain interval is either limited to this value or is the nearest to this
#'  value. A value <= .5 is useless.
#'@param spec.ui (default = .55). The specificity of the test scores within the
#'  uncertain interval is either limited to this value or is the nearest to this
#'  value. A value <= .5 is useless.
#'@param intersection (default = NULL) When NULL, the intersection is calculated
#'  with \code{get.intersection}, which uses the kernel density method to obtain
#'  the intersection. When another value is assigned to this parameter, this
#'  value is used instead.
#'@param return.first (default = TRUE) Return only the widest possible interval,
#'  given the restrictions. When FALSE all calculated intervals with their
#'  sensitivity and specificity are returned. NOTE: This function does not
#'  always find a suitable interval and can return a vector of NULL values.
#'@param select (default = 'nearest') If 'nearest', sensitivity and specificity
#'  of the uncertain interval are nearest sens.ui and spec.ui respectively. When
#'  'limited' the solutions have an uncertain interval with a sensitivity and
#'  specificity limited by sens.ui and spec.ui respectively.
#'
#'@details{ This function can be used for a test without a defined distribution
#'of the continuous test scores. The Uncertain Interval is generally defined as
#'an interval below and above the intersection, where the densities of the two
#'distributions of patients with and without the targeted condition are about
#'equal. This function uses for the definition of the uncertain interval a
#'sensitivity and specificity of the uncertain test scores below a desired value
#'(default .55).
#'
#'This essentially non-parametric function finds the best possible solution for
#'a sample. This function can be used for test with continuous scores or for
#'test with about twenty or more ordered test scores. The Uncertain Interval is
#'defined as an interval below and above the intersection, with a sensitivity
#'and specificity nearby or below a desired value (default .55).
#'
#'In its core, the \code{ui.nonpar} function is non-parametric, but it uses the
#'gaussian kernel for estimating the intersection between the two distributions.
#'Always check whether your results are within reason. If the results are
#'unsatisfactory, first check on the intersection. The \code{density} function
#'allows for other approximations than gaussian. Another estimate can be
#'obtained by using a more suitable kernel in the \code{density} function. The
#'parameter \code{intersection} can be used to assign the new estimate to the
#'\code{uncertain.interval} method.
#'
#'Furthermore, only a single intersection is assumed (or a second intersection
#'where the overlap is negligible). It should be noted that in most cases, a
#'test with more than one intersection with non-negligible overlap is
#'problematic and difficult to apply.
#'
#'The Uncertain interval method is developed for continuous distributions,
#'although it can be applied to ordered tests with distinguishable
#'distributions. When a test is used with less than 20 discernible values, a
#'warning is issued. The method may work satisfactorily, but results should
#'always be checked carefully.
#'
#'In general, when estimating decision thresholds, a sample of sufficient size
#'should be used. It is recommended to use at least a sample of 100 patients
#'with the targeted condition, and a 'healthy' sample (without the targeted
#'condition) of the same size or larger.
#'
#'The Uncertain interval method is not always capable to deliver results,
#'especially when select == 'limited'. Clearly, when there is no overlap between
#'the two distributions, there cannot be an uncertain interval. A very small
#'interval of overlap can also limit the possibilities to find a solution. When
#'there is no solution found, a vector of NA values is returned. }
#'@return {A \code{data.frame} of \describe{ \item{cp.l}{ Lower bound of the
#'  Uncertain interval.} \item{cp.h}{ Upper bound of the Uncertain interval.}
#'  \item{FN}{ Count of false negatives within the Uncertain interval.}
#'  \item{TP}{ Count of true positives within the Uncertain interval.}
#'  \item{TN}{ Count of true negatives within the Uncertain interval.}
#'  \item{FP}{ Count of false positives within the Uncertain interval.}
#'  \item{sensitivity}{ Sensitivity of the test scores within the Uncertain
#'  interval.} \item{specificity}{ Specificity of the test scores within the
#'  Uncertain interval.} } Only a single row is returned when parameter
#'  \code{return.first} = TRUE (default).}
#'@importFrom reshape2 melt
#'@references Landsheer, J. A. (2016). Interval of Uncertainty: An Alternative
#'  Approach for the Determination of Decision Thresholds, with an Illustrative
#'  Application for the Prediction of Prostate Cancer. PloS One, 11(11),
#'  e0166007.
#'
#'  Landsheer, J. A. (2018). The Clinical Relevance of Methods for Handling
#'  Inconclusive Medical Test Results: Quantification of Uncertainty in Medical
#'  Decision-Making and Screening. Diagnostics, 8(2), 32.
#'  https://doi.org/10.3390/diagnostics8020032
#'
#' @examples
#' # A simple test model
#' set.seed(1)
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,1))
#' ui.nonpar(ref, test, select='limited')
#'
#' ref = c(rep(0,20), rep(1,20))
#' test= c(rnorm(20), rnorm(20, mean=1))
#' ui.nonpar(ref, test)
#'
#'@export
ui.nonpar <-
  function(ref,
           test,
           sens.ui = .55,
           spec.ui = .55,
           intersection = NULL,
           return.first = T,
           select = c('nearest', 'limited')) {

    find.closest <- function(M, crit){
      mindiff=min(abs(M-crit))
      which((M == crit+mindiff) | (M == crit-mindiff), arr.ind=T)
    }

    bootstrap=0 # experimental
    select = match.arg(select)
    df = check.data(ref, test, model='kernel')

    if (bootstrap > 0) {
      # o4 = uncertain.interval(df$ref, df$test, sens.ui, spec.ui, intersection, return.first=T,
      #                         select, bootstrap=0)
      o4 = ui.nonpar(df$ref, df$test, sens.ui, spec.ui, intersection, return.first=T,
                              select)
      n = nrow(df) # n=10
      for (i in 1:bootstrap){
        sa = sample(n, replace=T)
        # o5 = uncertain.interval(df$ref[sa], df$test[sa], sens.ui, spec.ui, intersection, return.first=T,
        #                         select, bootstrap=0)
        o5 = ui.nonpar(df$ref[sa], df$test[sa], sens.ui, spec.ui, intersection, return.first=T,
                                select)
        o4  = rbind(o4,o5)
      }
      if (return.first){
        return(list(sample.est=o4[1,], boostrap.est=colMeans(o4[-1,])))
      } else {
        return(list(sample.est=o4[1,], boostrap.est=o4[-1,]))
      }
    }
    if (sens.ui < .5) stop('Value < .5 invalid for sens.ui')
    if (spec.ui < .5) stop('Value < .5 invalid for spec.ui')
    # if (sens.ui > .6) warning('Value > .6 not recommended for sens.ui')
    # if (spec.ui > .6) warning('Value > .6 not recommended for spec.ui')

    # only one relevant intersection assumed!
    # linear tests are used for determination of point of intersection
    # linear test is assumed to have a normal distribution
    if (is.null(intersection)) {
      intersection = get.intersection(df$ref, df$test, model='kernel')
      if (length(intersection) > 1) {
        intersection = utils::tail(intersection, n = 1)
        warning('More than one point of intersection. Highest density used. \n')
      } else intersection=intersection[1] # other values are ignored
    }

    tt <- table(df$test, df$ref) # sort test values; colSums(tt)

    pred = sort(unique(df$test))

    wi0 = rev(which(pred < intersection)) # count backwards from intersection
    wi1 = which(pred >= intersection)    # count forwards from intersection

    o0 = matrix(NA, ncol = 3, nrow = length(wi0))
    colnames(o0) = c('cp.l', 'TN', 'FN')
    temp = cumsum(tt[wi0, 1]) # colSums(tt[wi0[1:47],])
    o0[, 'cp.l'] = pred[wi0]  # as.numeric(names(temp))
    o0[, 'TN'] = temp
    o0[, 'FN'] = cumsum(tt[wi0, 2]) # head(o0); sprintf('%.20f', ua[1]) # include threshold
    # sprintf('%.20f', ua[1]); sprintf('%.20f', pred[wi0[47]])
    # ua[1] >=pred[wi0[47]]; which(pred >= ua[1]); which(abs(test-ua[1]) <= 1e-5)
    # pred[289]>=ua[1] ; test[237]>=ua[1];
    # sprintf('%.20f', pred[289]); sprintf('%.20f', test[237]); sprintf('%.20f', ua[1])

    o1 = matrix(NA, ncol = 3, nrow = length(wi1))
    colnames(o1) = c('cp.h', 'FP', 'TP')
    temp = cumsum(tt[wi1, 2])
    o1[, 'cp.h'] = pred[wi1] # as.numeric(names(temp))
    o1[, 'TP'] = temp
    o1[, 'FP'] = cumsum(tt[wi1, 1])  #head(o1)

    # find rows where TN/(TN+FP) <= sens.ui
    value = sens.ui / (1 - sens.ui)
    # r = o0[2, "TN"]
    if (select == 'limited') {
      res0 = lapply(o0[, "TN"], function(r) {
        a = which(r <= o1[, "FP"] * value)
        ifelse(length(a) == 0, return(NA), return(a) )})
      } else {
      res0 = lapply(o0[, "TN"], function(r) {
        a = find.closest(o1[, "FP"] * value, r)
        ifelse(length(a) == 0, return(NA), return(a) )} )
      }

    # create matrix from list
    m = t(sapply(res0, '[', 1:max(sapply(res0, length)))) # dim(m)
    # find rows with at least one valid finding
    rcp.l = which(rowSums(is.na(m)) != ncol(m))
    cp.l = o0[rcp.l, 'cp.l'] # candidates cp.l #cp.l=o0[m, 'cp.l'] # res0[[1]]

    # find candidates for cp.h
    cp.h = matrix(o1[m[rcp.l, ], 'cp.h'], nrow = length(cp.l)) # length(o1[m,'cp.h']) # length(cp.h)
    #  df=data.frame(cp.l, cp.h)
    df = data.frame(cbind(cp.l, cp.h))
    df.l = melt(df, id = 'cp.l')
    df.l = stats::na.omit(df.l)
    if (ncol(df.l) != 3)
      { oo1 =
      data.frame(
        'cp.l' = NA,
        'cp.h' = NA,
        'FN' = NA,
        'TP' = NA,
        'TN' = NA,
        'FP' = NA
      ) } else {
      colnames(df.l) <- c('cp.l', 'variable', 'cp.h')
      df.l = df.l[order(df.l$cp.l), c('cp.l', 'cp.h')]

      m.cp.l = match(df.l$cp.l, o0[, 'cp.l'])
      m.cp.h = match(df.l$cp.h, o1[, 'cp.h'])
      oo1 = data.frame(
        'cp.l' = df.l$cp.l,
        'cp.h' = df.l$cp.h,
        'FN' = o0[m.cp.l, 'FN'],
        'TP' = o1[m.cp.h, 'TP'],
        'TN' = o0[m.cp.l, 'TN'],
        'FP' = o1[m.cp.h , 'FP']
      )
      oo1$sensitivity = oo1$TP / (oo1$TP + oo1$FN)
      oo1$specificity = oo1$TN / (oo1$FP + oo1$TN) # head(oo1)
      if (select == 'limited') {
        oo1 = oo1[oo1$specificity <= spec.ui &
                    oo1$sensitivity <= sens.ui &
                    stats::complete.cases(oo1), ]
      } else {
        oo1 = oo1[stats::complete.cases(oo1), ]
      }
      oo1 = oo1[!duplicated(oo1[, c('FN', 'TP', 'TN', 'FP')]), ] # nrow(o1)
    }

    # find rows where TP/(TP+FN) <= spec.ui
    if (select == 'limited') {
      res1 = lapply(o1[, "TP"], function(r) {
        a = which(r <= o0[, "FN"] * spec.ui / (1 - spec.ui))
        ifelse(length(a) == 0, return(NA), return(a)) })
    } else {
      res1 = lapply(o1[, "TP"], function(r) {
        a = find.closest(o0[, "FN"] * value, r)
        ifelse(length(a) == 0, return(NA), return(a)) })
      }

    # create matrix from list
    m = t(sapply(res1, '[', 1:max(sapply(res1, length))))
    rcp.h = which(rowSums(is.na(m)) != ncol(m))
    cp.h = o1[rcp.h, 'cp.h'] # candidates cp.h

    # find candidates for cp.l
    cp.l = matrix(o0[m[rcp.h, ], 'cp.l'], nrow = length(cp.h))
    df = data.frame(cbind(cp.l, cp.h))
    df.h = melt(df, id = 'cp.h')
    df.h = stats::na.omit(df.h)
    if (ncol(df.h) != 3) {
      oo2 =
      data.frame(
        'cp.l' = NA,
        'cp.h' = NA,
        'FN' = NA,
        'TP' = NA,
        'TN' = NA,
        'FP' = NA
      )
    } else {
      colnames(df.h) <-
        c('cp.h', 'variable', 'cp.l') # error! Error in `colnames<-`(`*tmp*`, value = c("cp.h", "variable", "cp.l")) :
      # 'names' attribute [3] must be the same length as the vector [1]
      #df.h=df.h[df.h$cp.l > .25,]
      df.h = df.h[order(df.h$cp.h), c('cp.l', 'cp.h')]
      # colnames(df.h) <- c('cp.l', 'cp.h')

      m.cp.l = match(df.h$cp.l, o0[, 'cp.l'])
      m.cp.h = match(df.h$cp.h, o1[, 'cp.h'])
      oo2 = data.frame(
        'cp.l' = df.h$cp.l,
        'cp.h' = df.h$cp.h,
        'FN' = o0[m.cp.l, 'FN'],
        'TP' = o1[m.cp.h, 'TP'],
        'TN' = o0[m.cp.l, 'TN'],
        'FP' = o1[m.cp.h , 'FP']
      )
      oo2$sensitivity = oo2$TP / (oo2$TP +
                                    oo2$FN)
      oo2$specificity = oo2$TN / (oo2$FP +
                                    oo2$TN) # head(oo2)
      if (select == 'limited') {
        oo2 = oo2[oo2$specificity <= spec.ui &
                    oo2$sensitivity <= sens.ui & stats::complete.cases(oo2),]
      } else {
        oo2 = oo2[stats::complete.cases(oo2),]
      }
      oo2 = oo2[!duplicated(oo2[, c('FN', 'TP', 'TN', 'FP')]),] # nrow(o2)
    }
    o3 = rbind(oo1, oo2) # nrow(o3) # oo1=data.frame()
    o3 = o3[stats::complete.cases(o3), ]
    o3 = o3[order(-o3$cp.h, o3$cp.l),]
    o3 = o3[!duplicated(o3[, c('cp.l', 'cp.h')]), ]
    o3 = o3[!duplicated(o3[, c('FN', 'TP', 'TN', 'FP')]), ] # head(o3,10) # o=t(do.call(rbind, o3))
    if (select == 'nearest') {
      o3 = o3[find.closest(abs(o3[, 'sensitivity'] - sens.ui) +
                           abs(o3[, 'specificity'] - spec.ui), 0),]
    }
    if (return.first |
        nrow(o3) == 0)
      return(c(unlist(o3[1, ]), intersection=intersection))
    else
      return(t(do.call(rbind, o3)))
  }


