library(shiny)
library(UncertainInterval)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  # titlePanel("Medical Decision Methods: Uncertain Interval Determination for Binormal Distributed Test Scores"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    
    sidebarPanel(width = 6,
                 
                 plotOutput("distPlot"),        
                 
                 fluidRow(
                   p("The point of intersection is equal to the optimal dichotomous threshold at which the sum of Sensitivity and Specificity (Se + Sp) is maximized,
               and the sum of errors (FNR + FPR) is minimized. You manipulate the percentages of errors with the sliders."),
                   
                   column(6,
                          sliderInput("FNR",
                                      "False Negative Rate: Percentage of true patients with test scores below the intersection (1 - Se):",
                                      min = 1,
                                      max = 50,
                                      value = 15,
                                      #width='50%'
                          ),
                          textOutput("overlap"),
                          checkboxInput("combineSliders", 
                                        "Combine non-patients slider with patient slider", TRUE),
                   ), 
                   column(6,
                          sliderInput("FPR",
                                      "False Positive Rate: Percentage of true non-patients with test scores above the intersection (1 - Sp):",
                                      min = 1,
                                      max = 50,
                                      value = 15,
                                      #width = '50%'
                          ),
                          radioButtons("acc", label = h5("Allowed uncertainty (Se and Sp of Uncertain Interval):"),
                                       choices = list("Se.UI = Sp.UI = .55" = 1, "Se.UI = Sp.UI = .60" = 2),
                                       inline=F,
                                       selected = 1)
                   )
                 )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(width = 6,
              h3("Medical Decision Methods: Uncertain Interval Determination for Binormal Distributed Test Scores"),
              tabsetPanel(
                tabPanel("Introduction", 
                         h2('Determination of an Uncertain Interval'),
                         p("This Shiny App is intended to provide a hands-on demonstration of the",
                           span("binormal", style="color:blue"),  "Uncertain Interval Method for evaluating medical tests. 
             The Uncertain Interval method defines a range of test scores which 
             insufficiently distinguishes patients with and without the targeted disease. 
             Test scores outside this range 
             are used for a positive or negative classification, but for patients 
             with uncertain test outcomes, additional testing or waiting 
             on further developments (active  surveillance or watchful waiting) 
             is expected to be more beneficial.
             This Shiny App shows advantages and disadvantages of this evaluation 
             method of tests. It is intended to help you to decide for or against 
             the use of this method for evaluating your medical test."),
                         tags$ul(
                           tags$li("Hands-on. Look at the Hands-on tab for some examples and their discussion."),
                           tags$li("Control panel. Use this tab for explanation of the use of the controls."),
                           tags$li("Background. Explanation of the method."),
                           tags$li("Conclusion and References. My conclusions (but please draw your own) and References."),
                         ), 
                ), 
                
                tabPanel("Hands-on", 
                         h2("Hands-on demonstration of the bi-normal Uncertain Interval method"),
                         p("Set the uncertainty criterion to .55 and check the checkbox to 
             combine the two sliders. With the sliders you manipulate the overlap 
                           between the two distributions and the strength of the 
                           simulated test. A lower overlap directly defines a stronger test. For convenience, the",
                           span("AUC statistic", style="color:blue"), "is also presented, which is 
             a commonly used estimate of the strength of the test (using all 
                           test scores)."),
                         tags$ul(
                           tags$li("1. Make the test stronger by reducing the overlap between 
               the two distributions of the simulated test scores to two times 10%. The results are shown in the plot: the two 
               dashed vertical black lines that show the Uncertain Interval, the range of test 
               scores that hardly differentiate the two groups."),
                           p(span("UI.size", style="color:blue"), 'shows the total percentage 
                      of tested persons (true patients plus true non-patients) 
                      that have a test score that falls within the Uncertain Interval.'),
                           p(span("UI.err", style="color:blue"), 'shows the total percentage 
                      of classification errors avoided by recognizing the uncertain status in 
                      comparison to dichotomous optimal classification. Optimal dichotomous classification uses the 
                               optimal Youden threshold (that is the point of intersection). Please 
                               note that the percentage of errors avoided is considerably 
                               higher than the percentage of persons within the UI. (When the dischotomous optimal threshold is applied,
                               Se + Sp is maximized and the sum of errors is minimized. 
                               Using another dichotomous cutpoint would increase the number of errors.)'),
                           p('The values', 
                             span('Se.opt and Sp.opt', style="color:blue"), 
                             'show the optimal sensitivity and the specificity when the test 
                               scores are dichotomized using the intersection.',
                             span('Se.MCI and Sp.MCI', style="color:blue"), 
                             'show the sensitivity 
                               and specificity of the test scores outside the Uncertain Interval, 
                               when used for positive and negative classifications. 
                           The latter are always better and only equal in the case of 100% overlap.'),
                           p("A test with only 20% 
                       overlap is a relatively strong test, but it still has
                       uncertain test scores that are difficult to interpret:
                       test scores that have about equal probabilities to 
                       belong to either the group of patients or the group of non-patients."), 
                           tags$li("2. Make the test even stronger by reducing the percentages further. 
                       The Uncertain Interval becomes smaller when the test is stronger."), 
                           tags$li("3. Move the sliders to values larger than 10%. Observe that 
                      the Uncertain Interval grows when the test has more overlap 
                      and becomes weaker. Please note that the uncertain interval 
                      grows larger and larger when the test has too much 
                      overlap and is therefore of low quality. Using a weak test with a large 
                      uncertain interval is not very useful."),
                           tags$li("4. Move the sliders to 46. Please note that the reached 
                       values for Se.UI and Sp.UI (explained below) are now lower than the set 
                       value of .55. The values shown for ", 
                                   span('Sp.UI and Se.UI are the 
             values actually reached', style="color:blue"), 
                                   "for the test scores within the Uncertain Interval, 
             not the allowed uncertainty that has been set. Also note that the simulated test shows almost 
                       complete overlap for the scores of patients and 
                       non-patients. This simulated test is weak and in fact useless.")
                         ),
                         p("Please note that the values for Se.MCI and Sp.MCI are always larger than 
             Se.opt and Sp.opt. The meaning of this is that the 
                           Sensitivy and Specificity improves for the patients with a test score outside the Uncertain Interval,
                           that is, the patients who actually receive a positive or negative classification."),
                         p("Also note that the uncertain interval is always centered around 
                         the point of intersection, where the two distributions 
             have equal densities. The point of intersection is also the optimal 
             dichotomous cutpoint where the sum of Se + Sp is maximized and therefore 
             the sum of dischotomous classification errors (1 - Se + 1 - Sp) would be minimized. "),
                         p('In the examples above, both distributions 
           keep a standard deviation of 1. Remove the check from the checkbox to disconnect the two sliders. 
             Set the percentage of non-patients to 6 and the percentage of 
             patients to 26. The obtained distributions are now N(0, 1) and 
             N(3.31, 2.72). Observe that the Uncertain Interval remains centered around the 
             point of intersection.'),
                         tags$ul(
                           tags$li("5. When the test scores have different standard deviations 
               for patients and non-patients, the Uncertain Interval still indicates the test scores 
             that insufficiently distinguish the two groups. Earlier, this has been shown in simulation 
             studies as well with real data (Landsheer, 2016, 2018, 2020). Please 
             note the possibility of a secondary point of intersection when the 
             standard deviations for the two groups deviate. 
             For a test this is always undesirable, as the test 
             scores become inconsistent in their meaning for classification."),
                           p('Please note that the values shown for ', 
                             span('Sp.UI and Se.UI are the 
             values actually reached', style="color:blue"), 
                             'for the test scores within the Uncertain Interval, 
             not the allowed uncertainty that has been set.', 
                             'When the test score variance differ for both groups, the sum of Se.MCI and Sp.MCI is always greater than 
             the sum of Se.opt and Sp.opt; for weak tests they can be equal. N.B. Sp.MCI and Se.MCI are calculated in the same 
             way as Sp.VR and Se.VR for for the Virtual range of TG-ROC.'))
                ), 
                
                tabPanel("Control panel",           
                         h2('Using the dash board'),
                         p('The grey panel offers a dash board where the user can simulate many 
             different tests with bi-normal distributions. The tests differ in quality by 
             their overlap of the scores for the two groups. The overlap is 
             chosen with two sliders: The 
             upper slider sets the percentage of patients with test scores 
             above the point of intersection (blue dotted line in the left graph), 
             indicating the amount of False Negatives when the optimal dichotoumous 
             cut-point would be used. The 
             lower slider sets the percentage of patients which are known to not 
             have the disease and who have test scores below the point of intersection, 
             indicating the amount of False Negatives for the 
             dichotomous optimal cutpoint.'),
                         p('The true presence or absence of the 
             disease is known and is determined with superior means, called a "gold standard".'),
                         p('A checkbox allows the combinations of the two sliders. When the two sliders are combined, 
             the variance remains equal for 
             the two distributions. Unchecking makes it possible to use the two 
             sliders seperately, allowing the two distributions to have different variance.'),
                         p('In this application to demonstrate the working of the binormal 
             uncertain interval method, the criterion for uncertainty of the test scores 
             can be set to .55 or .60. The Uncertain Interval is always placed 
             around the point of intersection, which is also the optimal dichotomous 
             threshold where the Youden Index (Youden, 1950) or sum of 
             Sensitivity and Specificity is maximized.  The iterative method 
             defines in each iteration two borders where the Sensitivity and 
             Specificity of the test scores within these borders', 
                           span("(Sp.UI and Se.UI)", style="color:blue"), 'are 
             equal or lower than the citerion that is chosen (for instance, .55). The point of 
             intersection is used as the dichotomous threshold within this 
             Uncertain Interval. The iterative uncertain interval method stops when the criterion is reached. 
             The borders of the Uncertain Interval are therefore 
             always to the left and right of the point of intersection.'),
                         p('Consequently, the Sensitivity and Specificity of the Uncertain Interval are 
               directly connected to the optimal Sensitivity and Specificity:'), 
                         p(span("Se.opt = p * Se.UI + (1 - p) * Se.MCI", style="color:blue"),  ', and'),
                         p(span("Sp.opt = q * Sp.UI + (1 - q) * Sp.MCI", style="color:blue"),
                           ', where p and q are respectively the proportion of true patients and 
             non-patients within the Uncertain Interval. The values', 
                           span('Se.opt and Sp.opt', style="color:blue"), 
                           'show the optimal sensitivity and the specificity when the test 
             scores are dichotomized using the intersection.',
                           span('Se.MCI and Sp.MCI', style="color:blue"), 
                           'show the sensitivity 
             and specificity of the test scores outside the Uncertain Interval, 
             when these test scores are used for positive and negative classifications. '),
                         p("The total overlap between the two distributions defines 
                         directly the 
             strength of the simulated test. The overlap is here defined as the sum of the 
             two percentages set by the sliders. In this way, 
             a large amount of bi-normal tests of varying strenths can be simulated. 
             The strength of the test is directly determined by the overlap of 
             the distributions of test scores: a test is stronger when the 
             overlap is smaller. "),
                ),
                
                tabPanel("Background", 
                         h2("Background Information"),
                         p('A test, intended for comfirming the presence or absence of a 
             disease, can be evaluated using two groups: a group of true patients 
             who have the disease and a group of patients who truly do not have 
             the disease (shortly called non-patients). The plot shows 
             the two simulated densities of the the test scores of the two groups. 
             Higher test scores indicate the presence of the disease. The ', 
                           span('ui.binormal method from the R-package UncertainInterval', style="color:blue"), 
                           'determines iteratively the borders of the 
             Uncertain Interval (UI) for these test scores. It is assumed that the test scores are
             normally distributed. When this assumption cannot be made, other 
             methods from the UncertainInterval package may be more suitable. 
             This UI indcates a range of test scores that cannot sufficiently 
             distinguish the two groups of patients. '),
                         
                         h2('Test score Densities'),
                         p("The distributions shown in the plot show the densities 
             of the obtained simulated test scores, which are assumed to be normal 
             distributed. The test scores of the 
             non-patients are standard normal distributed, with mean of 0 
             and standard deviation of 1: N(0, 1). The distribution of the true 
             patients can vary widely. For a given test score, the difference of the two densities 
             indicates from which of the two groups of patients
             the test score is most likely drawn. When the sliders are combined, 
             both simualted distributions have a standard deviation of 1 and only the means 
             differs. The application starts with a distribution of test scores 
             with a mean of 2.07 and a standard deviation of 1 (N(2.07, 1))."),
                         
                         h2('Using Densities'),
                         p('In the Uncertain Interval method, densities or relative frequencies are 
                used for the determination whether a testcore (or a small range of test 
                scores) is sampled from the group of patients or from the 
                group of non-patients. Uncertain Interval method are intended for classification 
                of a patient, of whom it is suspected that he or she may 
                have the targeted disease. The pretest 
                probability is set to .50. In other words, before testing  it is considered 
                equally likely that the patient may or may not have the 
                disease. The method is therefore not fit for estimation of the probability
                that a person drawn from the general population has a particular 
                disease. For most diseases the prevalence of the disease 
                in the general propulation is much smaller than .5.'),
                         
                         h2("Trichotomization versus dichotomization"),
                         p("The classic techniques for evaluating tests for medical decision 
             make use of dichotomization of the test scores. Implicitely, all test scores are 
             considered as equally valid, and all test scores are used for a positive or 
             negative classification. The R package UncertainInterval covers 
             several trichotomization methods, that 
             try to identify test scores that are insufficiently valid and are 
             better not used for classification. The idea behind these trichotomization methods 
             is that patients are better served when these invalid test scores 
             are not used for a positive or negative classification. These 
             patients are better off with retesting, if 
             possible with a better test, or awaiting further developments. In 
             some medical fields the uncertainty of a diagnostic outcome is followed 
             up by techniques such as watchful waiting and active surveillance. 
             These more cautious approaches are especially relevant when possible 
             treatments can have serious side effects."),
                ),
                
                tabPanel("Conclusion and References", 
                         h2('In conclusion'),
                         p('The largest advantage of using the Uncertain Interval methods is 
           the identification of patients with uncertain, inconclusive test scores. 
           This allows for following a more cautious path for these patients. 
           Furthermore, the simulations show that the Uncertain Interval method for 
             bi-normal data almost always offer a better positive or negative 
             classification of 
             patients, with a higher sensitivity, respectively specificity,  
             when solely the test scores outside the Uncertain Interval are 
             used. When compared to using all test scores using the optimal 
             dichotomous cut-point (intersection or Youden threshold). In 
             Landsheer (2016, 2018, 2020) it has been shown that the Uncertain Interval 
             contain a relative large part of all classification errors and removing 
             these test scores from classification prevents a considerable amount 
             of errors when compared with  applying the dichotomous optimal  cut-point (Youden, 1950). 
             Further it may be noted that the Uncertain Interval 
             methods have none of the inconsistencies that plague the TG-ROC 
             method. '),
                         h2('References'),
                         p('Landsheer, J. A. (2016). Interval of uncertainty: an alternative 
             approach for the determination of decision thresholds, with an 
             illustrative application for the prediction of prostate cancer. 
             PloS one, 11(11), e0166007.',
                           tags$a(href="https://doi.org/10.1371/journal.pone.0166007", 
                                  "https://doi.org/10.1371/journal.pone.0166007")),
                         p('Landsheer, J. A. (2018). The clinical relevance of methods for 
             handling inconclusive medical test results: quantification of 
             uncertainty in medical decision-making and screening. 
             Diagnostics, 8(2), 32.', 
                           tags$a(href="https://doi.org/10.3390/diagnostics8020032", 
                                  "https://doi.org/10.3390/diagnostics8020032")),
                         p('Landsheer, J. A. (2020). Impact of the Prevalence of Cognitive 
             Impairment on the Accuracy of the Montreal Cognitive Assessment: 
             The Advantage of Using two MoCA Thresholds to Identify Error-prone 
             Test Scores. Alzheimer Disease & Associated Disorders.',
                           tags$a(href="https://dx.doi.org/10.1097/WAD.0000000000000365", 
                                  "https://dx.doi.org/10.1097/WAD.0000000000000365")),
                         p('Landsheer J. A. (2020). CRAN - Package UncertainInterval 0.6.0.',
                           tags$a(href=" https://CRAN.R-project.org/package=UncertainInterval", 
                                  " https://CRAN.R-project.org/package=UncertainInterval")),
                         p('Youden, W. J. (1950). Index for rating diagnostic tests. 
           Cancer, 3(1), 32–35.',
                           tags$a(href="https://doi.org/10.1002/1097-0142(1950)3:1<32::AID-CNCR2820030106>3.0.CO;2-3", 
                                  "https://doi.org/10.1002/1097-0142(1950)3:1<32::AID-CNCR2820030106>3.0.CO;2-3")),
                         
                )
              )
    )
  ))
