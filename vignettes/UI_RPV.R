## ---- include=FALSE------------------------------------------------------
options(tinytex.verbose = TRUE)

## ------------------------------------------------------------------------
library(UncertainInterval)
data('synthdata_NACC')
head(synthdata_NACC) 
nrow(synthdata_NACC) 

## ------------------------------------------------------------------------
m1 = synthdata_NACC[!is.na(synthdata_NACC$MOCATOTS.1) & !is.na(synthdata_NACC$ref.1), ]


## ------------------------------------------------------------------------
addmargins(table(m1$ref.1, m1$MOCATOTS.1, useNA = 'always'))
barplotMD(m1$ref.1, m1$MOCATOTS.1)

## ------------------------------------------------------------------------
quality.threshold(m1$ref.1, -m1$MOCATOTS.1, threshold = -25)

## ------------------------------------------------------------------------
get.intersection(m1$ref.1, -m1$MOCATOTS.1)

## ------------------------------------------------------------------------
quality.threshold(m1$ref.1, -m1$MOCATOTS.1, threshold = -23)

## ------------------------------------------------------------------------
t = addmargins(table(m1$ref.1, m1$center, useNA = 'always')) 
t = rbind(t, t[2,]/(t[2,]+t[1,]))
to = t[,c(order(t[5,1:30]),31:32)]
rownames(to) = c('0', '1','<NA>','Sum', 'prev')
round(to, 3)
center = colnames(to)[1:30] # sorted on sample prevalence

## ------------------------------------------------------------------------
indm = matrix(NA, 30, 9)
yt = rep(NA, 30)
for (i in 1:30) {
  # i=1
  ref = m1[m1$center == center[i], ]$ref.1 
  unique(ref)
  test = -m1[m1$center == center[i], ]$MOCATOTS.1
  N0 = length(test[ref == 0])
  N1 = length(test[ref == 1])
  indm[i, ] = c(N0, N1,
                quality.threshold(ref, test, threshold = -23, 
                   model = 'ordinal')$indices[c(1, 4:9)])
}
colnames(indm) = c('n0', 'n1', 'prev', 'Sp', 'Se', 'NPV', 'PPV', 'SNPV', 'SPPV')
rownames(indm) = 1:30
round(indm, 3)


## ------------------------------------------------------------------------
round(cor(indm[,'prev'], indm[,c('NPV', 'PPV', 'Sp', 'Se', 'SNPV', 'SPPV')]), 2)

## ------------------------------------------------------------------------
which(indm[,'Se'] < .7)
which(indm[,'Sp'] < .7)

## ------------------------------------------------------------------------
which(indm[,'NPV'] >= .8)
which(indm[,'PPV'] >= .8)
which(indm[,'PPV'] >= .8 & indm[,'NPV'] >= .8)

## ------------------------------------------------------------------------
ddiff = (m1$vdate.2 - m1$vdate.1)
summary(ddiff)
library(psych)
ICC(na.omit(cbind(m1$MOCATOTS.1, m1$MOCATOTS.2)))

## ------------------------------------------------------------------------
timesel = (ddiff >= 335) & (ddiff <= 395)
ICC(na.omit(cbind(m1$MOCATOTS.1[timesel], m1$MOCATOTS.2[timesel])))
# ICC(na.omit(cbind(m1$MOCATOTS.1[timesel], m1$MOCATOTS.2[timesel])), lmer=FALSE)

## ------------------------------------------------------------------------
RPV(m1$ref.1, m1$MOCATOTS.1, reliability = .86, roll.length = 5)

## ------------------------------------------------------------------------
class23 = as.numeric(m1$MOCATOTS.1 <= 23)
err = m1$ref.1 != class23
sum(err[m1$MOCATOTS.1 >= 22 & m1$MOCATOTS.1 <= 25])/ sum(err)

## ------------------------------------------------------------------------
indm2 = matrix(NA, 30, 4); i=1
for (i in 1:30){
  ref = m1[m1$center==center[i],]$ref.1
  test = m1[m1$center==center[i],]$MOCATOTS.1 # reversed order
  res = RPV(ref, test, pretest.prob = .53, reliability=.86, roll.length=5,
            decision.odds = 2, preselected.thresholds = c(25,22))$result[4:6,]
  res2=c(t(res))[c(1,3,4,9)]
  indm2[i,] = as.numeric(c(sub("%","",res2[1]), sub("%","",res2[2]),
                          sub("%","",res2[3]),sub("%","",res2[4])))/100
}
indm2=cbind(to['prev',1:30], indm2)
colnames(indm2) = c('prev', 'NPV', 'PPV', 'Sp', 'Se')
# SPPV = Se / (Se + 1 – Sp) and SNPV = Sp / (Sp + 1 – Se)).
SNPV = indm2[,'Sp']/(indm2[,'Sp']+ 1 - indm2[,'Se'])
SPPV = indm2[,'Se']/(indm2[,'Se']+ 1 - indm2[,'Sp'])
rownames(indm2)= 1:30
indm2 = cbind(indm2, SNPV, SPPV)
round(indm2, 3)

## ------------------------------------------------------------------------
round(cor(indm2[,'prev'], 
          indm2[,c('NPV', 'PPV', 'Sp', 'Se', 'SNPV', 'SPPV')]), 2)

## ------------------------------------------------------------------------
which(indm2[,'NPV'] >= .8)
which(indm2[,'PPV'] >= .8)
which(indm2[,'PPV'] >= .8 & indm2[,'NPV'] >= .8)

