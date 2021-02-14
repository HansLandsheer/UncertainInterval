library(shiny)
# library(UncertainInterval)

# runExample()
server <- function(input, output, session) {
  
  output$distPlot <- renderPlot({
    
    # input=data.frame(
    # combineSliders = F,
    # FNR = 15,
    # FPR = 25,
    # pretestprob = 50,
    # criterion = 1)
    
    if (input$combineSliders) updateSliderInput(session, "FPR", value = input$FNR)
    
    pretestprob = prevalence = input$pretestprob 
    FPR = input$FPR / 100
    FNR = input$FNR / 100
    negpost.crit = input$negpost.crit
    pospost.crit = input$pospost.crit
    
    m0=0; sd0=1; 
    is = qnorm((1-FPR)) # x value for 15% FP = intersection
    Z = qnorm(FNR) # Z value = 15% FN; Z = (x-mean) / sd
    dens = dnorm(is)  * (1-prevalence) # density at the point of intersection
    sd1 = (prevalence/(dens*sqrt(2*pi)))*exp(-.5 * Z^2)
    m1 = is-Z*sd1 
    
    x <- seq(-4, 4+m1, length=1000)
    y0 <- dnorm(x, m0, sd0) * (1-prevalence)
    y1 = dnorm(x, m1, sd1) * prevalence
    
    par(mfrow=c(1,2))
    plot(x, y0, type="l", col='black', xlab='Test score', ylab='Density', 
         main='Bi-normal Grey Zone',
         sub = paste('N(', m0,', ', sd0,')','  N(', round(m1,2),', ', round(sd1,2),')', sep=''),
         ylim=c(0,.25), font.sub=3) 
    lines(x, y1, type="l", col='red')
    # position and density of point of intersection
    segments(is,0,is, dens, col='blue') 
    
    sp = pnorm(x, m0, sd0)
    se = 1 - pnorm(x, m1, sd1) # se and sp for all x-values
    
    preodds=pretestprob/(1-pretestprob); # pretest odds
    neg.lr=(1-se)/sp; # negative likelihood ratio
    pos.lr=se/(1-sp); # positive likelihood ratio
    # By Bayes' rule : posttest(posterior)odds=likelihoodratio * pretest(prior)odds
    negpostodds=neg.lr * preodds ; # negative post odds
    pospostodds=pos.lr * preodds ; # positive post odds
    negpostprob = negpostodds / (negpostodds+1); # negative post probability
    pospostprob=pospostodds / (pospostodds+1); # positive post probability
    
    # Find points where pospostprob is above pospost.crit.
    # https://stackoverflow.com/questions/20519431/finding-point-of-intersection-in-r
    above <- pospostprob > pospost.crit
    # Points always intersect when above=TRUE, then FALSE or reverse
    intersect.points.pos <- which(diff(above) != 0)
    # Find points where pospostprob is above pospost.crit.
    above <- negpostprob > negpost.crit
    # Points always intersect when above=TRUE, then FALSE or reverse
    intersect.points.neg <- which(diff(above) != 0)
    # determine points of intercept with criteria
    pos.border = x[intersect.points.pos[1]] # always use first value
    neg.border = x[intersect.points.neg[length(intersect.points.neg)]] # always use last value
    
    linetype.pos = c(1,2); if(length(intersect.points.pos) == 1) linetype.pos =1
    abline(v=x[intersect.points.pos], lty=linetype.pos, col='red')
    linetype.neg = c(1,2); if(length(intersect.points.neg) == 1) linetype.neg =1
    if(length(intersect.points.neg) == 1) {
      abline(v=x[intersect.points.neg], lty=1, col='black')
    } else {
      abline(v=x[intersect.points.neg], lty=c(2,1), col='black')
    }
    
    a = (m1-m0)/sd1
    b = sd0/sd1
    AUC = round(unname(pnorm(a/sqrt(1+b^2))), 4)
    positives.grey = prevalence * (pnorm(pos.border, m1, sd1) - pnorm(neg.border, m1, sd1)) # all grey scores positive
    negatives.grey = (1 - prevalence) * (pnorm(pos.border, m0, sd0) - pnorm(neg.border, m0, sd0)) # all grey scores negative
    # is = optimal cutpoint (intersection)
    FN.grey = pnorm(is, m1, sd1) - pnorm(neg.border, m1, sd1)
    FP.grey =  pnorm(pos.border, m0, sd0) - pnorm(is, m0, sd0)
    # err.grey = FP.grey + FN.grey
    legend('topleft', c(paste('AUC = ', AUC, sep=''),
      paste('Grey.size = ', round((positives.grey+negatives.grey)*100,0), '%', sep='')),
     # paste('Grey.err = ', round((err.grey/(FNR+FPR))*100, 0), '%', sep='')),
     # text.col= c('black','red', 'black', 'black')
           )
    legend('topright', 'intersection', lty=1, col=c('blue'), 
    #        text.col=c('red', 'black', 'blue', 'black')
    )
    
    Se.Grey = Sp.Grey = NA
    if (!length(pos.border) == 0 & !length(neg.border) == 0) {
      if (!is.na(pos.border) & !is.na(neg.border)) {
        if (pos.border > is & neg.border < is) {
          Se.Grey = (pnorm(pos.border, m1, sd1)-pnorm(is, m1, sd1)) / (pnorm(pos.border, m1, sd1) - pnorm(neg.border, m1, sd1))
          Sp.Grey = (pnorm(is, m0, sd0)- pnorm(neg.border, m0, sd0)) / (pnorm(pos.border, m0, sd0)- pnorm(neg.border, m0, sd0))
        }
      }
    }
     
    Se.opt = 1-pnorm(is, m1, sd1)
    Sp.opt = pnorm(is, m0, sd0)
    Se.nonGrey = (1-pnorm(pos.border, m1, sd1))/(pnorm(neg.border, m1, sd1)+ 1-pnorm(pos.border, m1, sd1))
    Sp.nonGrey = (pnorm(neg.border, m0, sd0))/(1-pnorm(pos.border, m0, sd0)+ pnorm(neg.border, m0, sd0))
    
    if (identical(Sp.Grey, numeric(0))) Sp.Grey = NA
    if (identical(Se.Grey, numeric(0))) Se.Grey = NA
    if (identical(Sp.nonGrey, numeric(0))) Sp.nonGrey = NA
    if (identical(Se.nonGrey, numeric(0))) Se.nonGrey = NA
    
    # legend('right', c(paste('Se.opt =', round(Se.opt,2)),
    #                   paste('Se.nonGrey =', round(Se.nonGrey,2))), text.col=c('red', 'red'))
    # legend('left', c(paste('Sp.opt =', round(Sp.opt,2)),
    #                  paste('Sp.nonGrey =', round(Sp.nonGrey,2))), text.col=c('black', 'black'))

    plot(x, pospostprob, type='l', col='red', ylim=c(0,1),
         main="Post probability distributions",
         xlab="Test score", ylab="Post Probability") 
    abline(h = pospost.crit, v = x[intersect.points.pos], col='red', lty=c(1,linetype.pos))
    lines(x, negpostprob, type='l')
    abline(h = negpost.crit, v = x[intersect.points.neg], col='black', lty=linetype.neg)
    
    legend("topleft", c('Pos. Post Prob.',
                        'Neg. Post Prob'), col=c('red', 'black'),lty=c(1,1))
    legend('bottomright', c("Neg. border Grey ", "Pos. border Grey"), lty=c(1, 1), 
                            col=c('black', 'red'))
    
    Se.table = data.frame(Sp=c(Sp.opt,Sp.Grey, Sp.nonGrey), Se=c(Se.opt,Se.Grey, Se.nonGrey), 
                    row.names = c("Optimal dichotomous results*", "Grey interval*", "Non-grey interval"))
    
    prevalence = pretestprob
    NPV.all = (1-prevalence)*pnorm(is, m0, sd0)/((1-prevalence)*pnorm(is, m0, sd0)+prevalence*pnorm(is, m1, sd1))
    PPV.all = prevalence*(1-pnorm(is, m1, sd1))/((1-prevalence)*(1-pnorm(is, m0, sd0))+prevalence*(1-pnorm(is, m1, sd1)))
    SNPV.all = pnorm(is, m0, sd0)/(pnorm(is, m0, sd0)+pnorm(is, m1, sd1))
    SPPV.all = (1-pnorm(is, m1, sd1))/(1-pnorm(is, m0, sd0)+(1-pnorm(is, m1, sd1)))
    negpost.all = 1-NPV.all
    pospost.all = PPV.all

    # it is assumed that the prevalence is equally distributed in the grey interval and in the non-grey intervals
    NPV.grey=PPV.grey=SNPV.grey=SPPV.grey=negpost.grey=pospost.grey=NA
    if ((length(pos.border > 0) & length(neg.border) > 0)) {
        if (!is.na(neg.border) & !is.na(pos.border)) { 
        positives.grey = pnorm(neg.border, m1, sd1) - pnorm(pos.border, m1, sd1) # all grey scores positive
        negatives.grey = pnorm(neg.border, m0, sd0) - pnorm(pos.border, m0, sd0) # all grey scores negative
        NPV.grey = (1-prevalence)*negatives.grey/((1-prevalence)*negatives.grey+prevalence*positives.grey)
        PPV.grey = (prevalence*(positives.grey))/((1-prevalence)*(negatives.grey)+prevalence*(positives.grey))
        SNPV.grey = negatives.grey/(negatives.grey+positives.grey)
        SPPV.grey = positives.grey/(negatives.grey+positives.grey)
        negpost.grey = 1-NPV.grey
        pospost.grey = PPV.grey
      }
    }
    # it is assumed that the prevalence is equally distributed in the grey interval and in the non-grey intervals
    # lng = lower non-grey interval
    NPV.lng=PPV.lng=SNPV.lng=SPPV.lng=negpost.lng=pospost.lng=NA
    if (length(neg.border > 0)) {
      if (!is.na(neg.border)) {
        positives.lng = pnorm(neg.border, m1, sd1) # all lng scores positive
        negatives.lng = pnorm(neg.border, m0, sd0) # all lng scores negative
        NPV.lng = (1-prevalence)*negatives.lng/((1-prevalence)*negatives.lng+prevalence*positives.lng)
        PPV.lng = (prevalence*(positives.lng))/((1-prevalence)*(negatives.lng)+prevalence*(positives.lng))
        SNPV.lng = negatives.lng/(negatives.lng+positives.lng)
        SPPV.lng = positives.lng/(negatives.lng+positives.lng)
        negpost.lng = 1-NPV.lng
        pospost.lng = PPV.lng 
      }
    }
    # ung = upper non-grey interval
    NPV.ung=PPV.ung=SNPV.ung=SPPV.ung=negpost.ung=pospost.ung=NA
    if (length(pos.border > 0)) {
      if (!is.na(pos.border)) {
        positives.ung = 1 - pnorm(pos.border, m1, sd1) # all ung scores positive
        negatives.ung = 1 - pnorm(pos.border, m0, sd0) # all ung scores negative
        NPV.ung = (1-prevalence)*negatives.ung/((1-prevalence)*negatives.ung+prevalence*positives.ung)
        PPV.ung = (prevalence*(positives.ung))/((1-prevalence)*(negatives.ung)+prevalence*(positives.ung))
        SNPV.ung = negatives.ung/(negatives.ung+positives.ung)
        SPPV.ung = positives.ung/(negatives.ung+positives.ung)
        negpost.ung = 1-NPV.ung
        pospost.ung = PPV.ung
      }
    }

        PV.table = data.frame(NPV=c(NPV.all,NPV.grey,NPV.lng,NPV.ung), 
                              PPV= c(PPV.all,PPV.grey,PPV.lng,PPV.ung), 
                              SNPV= c(SNPV.all,SNPV.grey,SNPV.lng,SNPV.ung), 
                              SPPV= c(SPPV.all,SPPV.grey,SPPV.lng,SPPV.ung), 
                              # "neg.post.pred.prob"= c(negpost.all,negpost.grey,negpost.lng,negpost.ung), 
                              "post.pred.prob"= c(pospost.all,pospost.grey,pospost.lng,pospost.ung),
                          row.names=c('All optimal dichotomous results*', 'Grey interval', 
                                      'Neg. non-grey interval (lower)', 'Pos. non-grey interval (upper)'))
    
    if (input$outputTable==1)  {
      ot = Se.table
      } else{
      ot = PV.table  
      }
    output$table = renderTable(ot, rownames=T)
    #  par(mfrow=c(1,1))
  })
  
  output$overlap = renderText(paste("Overlap Distributions =", 
      2*((1-input$pretestprob/100)*input$FNR+
           (input$pretestprob/100)*input$FPR), '%'))
  
}
