library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  # titlePanel("Medical Decision Methods: TG-ROC for Binormal Distributions"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      p("The point of intersection is equal to the optimal dichotomous threshold at which the sum of Sensitivity and Specificity (Se + Sp) is maximized,
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
      checkboxInput("combineSliders", "Combine non-patients slider with patient slider", TRUE),
      radioButtons("acc", label = h3("Accuracy level Se and Sp:"),
                   choices = list(".90" = 1, ".95" = 2), 
                   selected = 1)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot"),
      h1('Demonstration of Two-Graph Receiver Operating Characteristic (TG-ROC)'),
      
      h2("Hands-on Demonstration of TG-ROC for the determination of an Intermediate Range"),
      p('Set the accuracy level to .9 and check the checkbox to combine the two sliders.'),
      p("The two vertical dashed lines show the borders of the Intermediate 
             Range, associated with Se = .9 (red) and with Sp = .9 (black). 
             It is assumed that the scores within this Intermediate Range are inconclusive."),
      
      tabsetPanel(
        tabPanel("Introduction"),
        h2('Determination of an Intermediate Range')),
      tags$ul(
        tags$li("1. Make the test stronger by reducing the overlap between 
               the two distributions to two times 10%. The result is that the line for Se = .9 coincides with
                       the line for Sp = .9. In that case there is no 
                       Intermediate Range of inconclusive test csores, which is 
                       unexpected. A test with 20% 
                       overlap is relatively strong, but it still has
                       intermediate test scores that are difficult to interpret. That is, 
                       test scores that have about equal probabilities to 
                       be sampled from either the group of patients or the group of non-patients."), 
        tags$li("2. Make the test even stronger by reducing the percentages further. 
                       The Intermediate Range grows 
                       again, and this range grows larger as the test is stronger.
                       This is not what we want: a stronger test should have 
                       a smaller intermediate range of inconclusive results, not a larger one."), 
        tags$li("3. Move the sliders to values larger than 10%. Observe that 
                       the lines for Se and Sp change position. Now, the range of test 
                       scores that offer a Sensitivity of .9 INCLUDE the test 
                       scores of the intermediate range. Similarly, the test 
                       scores that offer a Specificity of .9 also INCLUDE the 
                       test scores of the Intermediate Range. This is confusing 
                       and can lead to over estimation of the Sensitivity and 
                       Specificity of the Valid Ranges of these weaker tests. In 
                       these cases, Se.VR and Sp.VR (explained below) are always lower than the chosen levels.")),
      p("There is yet another issue. In the examples above, both distributions 
           keep a standard deviation of 1. With equal variance,
             the Intermediate Range is centered around the point of intersection of the 
             two distributions. The point of intersection is where the two distributions 
             have equal densities. When the tests have different standard deviations, the Intermediate 
             Range moves to the left or the right of the point of intersection."),
      p('Remove the check from the checkbox to disconnect the two sliders. 
             Set the percentage of true patients to 26 and the percentage of 
             non-patients to 26. The obtained distributions are now N(0, 1) and 
             N(3.31, 2.72). The intermediate range is moved to the left and the 
             point of intersection falls outside the Intermediate Range.'),
      tags$ul(
        tags$li("4. When the point of intersection falls outside the Intermediate Range,
             some test scores within the intermediate range 
              have densities for the two distributions with a considerable difference. 
             In that situation, the Intermediate Range contains test scores 
             that are relatively easy to distinguish. This results in a worsening of the classification
             which is not desirable. Earlier, this has been demonstrated in a simulation 
             study (Landsheer, 2016). ")
      ),
      
      h2('Using the dash board'),
      p('The grey panel offers a dash board where the user can create many 
             different tests. The tests differ in their overlap of the scores 
             for the two groups. The overlap is chosen with two sliders: the 
             upper slider sets the percentage of patients with test scores 
             below the point of intersection (the intersection is the blue dotted 
             line in the left graph). The 
             lower slider sets the percentage of non-patients whith test scores above the intersection. The true presence or absence of the 
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
      p('When a test intended for comfirming the presence or absence of a 
             disease, the test is evaluated using two groups: a group of true patients 
             who have the disease and a group of patients who truly do not have 
             the disease (shortly called non-patients). The left plot above shows 
             the two densities of the two groups. The right plot shows the TG-ROC, 
             and shows the Sensitivity and Specificity when a test score is used as a 
             dichotomous cut-point for classifying all patients positive or negative 
             for the presence of the disease. The chosen accuracy level (.9 or .95) 
             is the desired level for Sensitivity and Specificity. Clearly, 
             Sensitivity increases while Specificity decreases and vice versa.'),
      
      h2('Two normal densities and TG-ROC'),
      p("The bi-normal distributions shown in the left plot show the densities 
             of the obtained simulated test scores. The test scores of the 
             non-patients are always standard normal distributed, with mean of 0 
             and standard deviation of 1: N(0, 1). The distribution of the true 
             patients can vary widely. The difference of the two densities 
             indicates for a given test score from which of the two groups of patients
             the test score is most likely drawn. When the sliders are combined, 
             both distributions have a standard deviation of 1 and only the means 
             differs. The application starts with both sliders set to 15%, which 
             results in a distribution of true patients' test scores 
             with a mean of 2.07 and a standard deviation of 1 (N(2.07, 1))."),
      p("TG-ROC shows the intersection of Sensitivity and Specificity with the chosen 
             accuracy level, in this case .9 or .95. The two dashed vertical lines 
             show the upper and lower border of the Intermediate range; they 
             are shown in both the left and right graph and are always the same in the two graphs.
             These two dashed lines represent the chosen accuravcy level 
             and indicate the two cut-off values, which are the
             'lower' and 'upper limits' of the Intermediate Range (IR). The scores 
             in this Intermediate Range are considered as inconclusive, 
             'non-positive, non-negative' test results. 
             According to Greiner (1995), 
             'Considering only results outside the IR ... the test's Se and Sp 
             would be 95 or 90%, respectively' (p.125). 
             Regretfully, this is not true.",
        span("Se.VR", style="color:blue"), "and", span("Sp.VR", style="color:blue"), 
        "show the realised sensitivity 
             and specificity of the test scores in the Valid Range, that includes 
             all test scores outside the Intermediate Range. These test scores in 
             the Valid Range are used 
             for a positive or negative classification. Although in some cases the values 
             for Se.VR and Sp.Vr are higher than the chosen accuracy level, but 
             in most cases they are lower. 
             "),
      
      h2("Trichotomization versus dichotomization"),
      p("The classic techniques for evaluating tests for medical decision 
             make use of dichotomization of the test scores. All test scores are 
             considered as equally useful for classification and all are used for a positive or 
             negative classification of each patient concerning the disease that 
             is the target of the test. TG-ROC is a trichotomization method, that 
             tries to identify test scores that are insufficiently valid and are 
             better not used for classification. Before classification, a range of the least 
             valid test scores is determined. The idea behind these trichotomization methods 
             is that patients are better served when these invalid test scores 
             are not used for classification. These patients are better off with retesting, if 
             possible with a better test, or awaiting further developments. In 
             some medical fields the uncertainty of a diagnostic outcome is followed 
             up by techniques such as watchful waiting and active surveillance. 
             These more cautious approaches are especially relevant when possible 
             treatments can have serious side effects."),
      
      h2('In conclusion'),
      p('The TG-ROC method for tricotomizations shows several inconsistencies. 
              TG-ROC does not offer the accuracies that are promised. For some 
             tests, an Intermediate Range with the chosen accuracy is not existent. 
             In other cases the Intermediate Range of supposedly inconclusive test 
             scores can be larger for stronger tests and smaller for 
             weaker tests. When the standard deviations of the results differ 
             for patients and non-patients, the Intermediate Range may exclude 
             test scores that offer a fairly good distinction between the two groups. 
             Use of TG-ROC can therefore not be recommended for identifying inconclusive test scores.'),
      h2('References'),
      p('Greiner, M., Sohr, D., & Göbel, P. (1995). A Modified ROC Analysis 
           for the Selection of Cut-Off Values and the Definition of Intermediate 
           Results of Serodiagnostic Tests. Journal of Immunological Methods, 
           185(1), 123–132.'),
      p('Landsheer, J. A. (2016). Interval of uncertainty: an alternative 
             approach for the determination of decision thresholds, with an 
             illustrative application for the prediction of prostate cancer. 
             PloS one, 11(11), e0166007.')
      
    )
  )
)