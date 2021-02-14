library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Medical Decision Methods: Receiver Operating Characteristics (ROC) for Binormal Distributions"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            p("The point of intersection is equal to the optimal dichtomous threshold at which the sum of Sensitivity and Specificity (Se + Sp) is maximized,
               and the sum of errors (FNR + FPR) is minimized. You manipulate the percentages of errors with the sliders."),
            sliderInput("FNR",
                        "False Negative Rate: Percentage of true patients with test scores below the intersection (1 - Se):",
                        min = 1,
                        max = 50,
                        value = 15),
            sliderInput("FPR",
                        "False Positive Rate: Percentage of true non-patients with test scores above the intersection (1 - Sp):",
                        min = 1,
                        max = 50,
                        value = 15),
            checkboxInput("combineSliders", "Combine patients slider with non-patient slider", TRUE),
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
           h1('Demonstration of the Receiver Operating Characteristic curve (ROC)'),

           h2("Hands-on Demonstration of ROC"),
           p("1. Check the checkbox to combine the two sliders. Move the upper 
           slider to the left or right, simulating a better or 
             worse test. ROC curves closer to the top-left corner indicate 
             better performance (AUC). A curve close to the diagonal indicates 
             a worthless test (AUC = .5). The ROC curves shows us the relative 
             trade-off between true positives (benefits) and false positives (costs). 
             The ROC curve does not show us the actual cut-off scores."),
           p("2. Uncheck the checkbox to release the lower slider. Move the 
             upper slider to 30, and simulate different tests with the lower slider. 
             The two tests have now different variances. Observe that the ROC 
             curve sometimes crosses the daigonal. In the left plot you can see 
             that this is caused by a secondary point of intersection. This is 
             undesirable, because the test scores around the second intersection 
             have a problematic and inconsistent interpretation: they no longer 
             allow us to say that only the higher scores indicate the targeted disease. "),
           p("3. Create different tests and observe the values of AUC. AUC is 
             abbreviation of 'Area under the Curve'. If you were wondering 'Which curve?', 
             well, it is the Area under the ROC curve, and AUROCC would be a better name. 
             As we do not want the curve 
             below the diagonal, AUC varies in practice between .5 
             (complete overlap) and 1.0 (no overlap)."),
 

           h2('Using the dash board'),
           p('The grey panel offers a dash board where the user can create many 
             different tests. The tests differ in their overlap of the scores 
             for the two groups. The overlap is chosen with two sliders: the 
             upper slider sets the percentage of true patients with test scores 
             below the point of intersection (that is the blue dotted line in the left graph). The 
             lower slider sets the percentage of patients which are known to not 
             have the disease that have test scores above the intersection. 
             The true presence or absence of the 
             disease is known and is determined with superior means, called a "gold standard".'),
           p('A checkbox allows the combinations of the two sliders. When the two sliders are combined, the variance remains equal for 
             the two distributions. Unchecking makes it possible to use the two 
             sliders seperately, allowing the two distributions to have different variance.'),
           p('The total overlap is here defined as the sum of these two percentages. In this way, 
             a large amount of tests of varying strenths can be simulated. 
             The strength of the test is directly determined by the overlap of 
             the distributions of test scores: a test is stronger when the 
             overlap is smaller. For convenience, the',
             span("AUC statistic", style="color:blue"), 'is presented, which is 
             also an estimate of the strength of the test.  '),
           
           h2("Background Information"),
           p("When a test intended for comfirming the presence or absence of a 
             disease, the test is evaluated using two groups: a group of true patients 
             who have the disease and a group of patients who truly do not have 
             the disease (shortly called non-patients). The true status of each 
             patient is determined with a 'gold standard'. The left plot above shows 
             the two densities of the two groups. The right plot shows the ROC, 
             and shows the trade-off between Sensitivity and Specificity. Clearly, 
             Sensitivity increases while Specificity decreases and vice versa."),

           h2('Two normal densities and ROC'),
           p("The bi-normal distributions shown in the left plot show the densities 
             of the obtained simulated test scores. The test scores of the 
             non-patients are always standard normal distributed, with mean of 0 
             and standard deviation of 1: N(0, 1). The distribution of the true 
             patients can vary widely. The difference of the two densities 
             indicates for a given test score from which of the two groups of patients
             the test score is most likely drawn. When the sliders are combined, 
             both distributions have a standard deviation of 1 and only the means 
             differs. The application starts with a distribution of test scores 
             with a mean of 2.07 and a standard deviation of 1 (N(2.07, 1))."),
           
           h2("Trichotomization versus dichotomization"),
           p("The classic techniques such as the ROC curve for evaluating tests for medical decision 
             make use of dichotomization of the test scores. All test scores are 
             considered as equally useful for a positive or 
             negative classification of each patient concerning the disease that 
             is the target of the test. Trichotomization methods say that some test 
             scores offer an insufficient distinction between the two groups of patients and 
             try to identify test scores that are insufficiently valid and are 
             better not used for classification. "),

           h2('In conclusion'),
           p('The ROC method is useful for showing the trade-off between 
           Sensitivity and Specificity. It is also useful to compare different 
           tests: a test is better when the curve is more drawn to the upperleft 
           corner. This is directly related to AUC: the Area Under the curve. 
           Furthermore, it allows us to identify short-comings of the test, as 
           it shows a crossing of the diagonal when the test has an undesirable 
           secondary point of intersection.'),
           p("The ROC curve shows us the trade-off between Sensitivity and 
             Specificity. These statistics are however group statistics. 
             Sensitivity gives us the percentage patients with a correct positive 
             classification and concerns the patients with test scores higher 
             than the cut-off score. Specificity gives us the percentage patients 
             with a correct negative classification and concerns the patients 
             with test scores higher than the cut-off score. Individual patients 
             do not have a true test score that is equal or higher than the cut-off 
             score, but their true score lies around the received test score. 
             For such a small interval of true scores around the received test 
             score, the left plot is easier 
             to interpret: simply look at the difference between the two density 
             plots. A large difference indicates good discrimination for these 
             test scores, a small difference indicates test scores that perhaps 
             are better not used for classification.
             "),
               
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    output$distPlot <- renderPlot({
        # generate pdf's based on input$bins from ui.R
        if (input$combineSliders) updateSliderInput(session, "FPR", value = input$FNR)
        FPR = input$FPR / 100
        FNR = input$FNR / 100
        acc = ifelse(input$acc==1, .9, .95)

        m0=0; sd0=1; 
        is = qnorm(1-FPR, 0, 1) # x value for 15% FP = intersection
        Z = qnorm(FNR) # Z value = 15% FN; Z = (x-mean) / sd
        dens = dnorm(is, 0, 1) # density at the point of intersection
        sd1 = (1/(dens*sqrt(2*pi)))*exp(-.5 * Z^2)
        m1 = is-Z*sd1 
        
        x <- seq(-4, 4+m1, length=200)
        y0 <- dnorm(x, m0, sd0)
        y1 = dnorm(x, m1, sd1)
        par(mfrow=c(1,2))
        plot(x, y0, type="l", col='black', xlab='Test score', ylab='Density', 
             main='Probability Density Functions', ylim=c(0,.5))
        lines(x, y1, type="l", col='red')
        lines(x=c(is,is), y=c(0,.5),col='blue', lty=3)
        a = (m1-m0)/sd1
        b = sd0/sd1
        AUC = round(unname(pnorm(a/sqrt(1+b^2))), 4)
        legend('topleft', c(paste('Non-Patients = N(', m0,', ', sd0,')', sep=''),
                            paste('True Patients = N(', round(m1,2),', ', round(sd1,2),')', sep='')),
                            text.col= c('black','red'))
        
        legend('topright', c('intersection',paste('AUC = ', AUC, sep='')), lty=c(3, 0),col=c('blue', 'black'))
        
        Se = 1-pnorm(x, m1, sd1)
        Sp = pnorm(x, 0, 1)
        plot(1-Sp, Se, type='l', col='red', xlab='1 - Specificity', 
             ylab='Sensitivity', main="ROC curve of 1-Specificity and Sensitivity") 
        lines(x=c(0,1), y=c(0,1),col='black') # diagonal
        par(mfrow=c(1,1))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
