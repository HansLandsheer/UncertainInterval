library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(sidebarLayout(
  sidebarPanel(
    sliderInput(
      "pretestprob",
      "Pre-test probability",
      min = 0.01,
      max = 1,
      value = .50,
      step=.01,
      width='100%'
    ),
    
    div(style = "float: left; width: 48%", # this number can be any percentage you need
        sliderInput(
          "FNR",
          "False Negative Rate: Percentage of true patients with test scores below the intersection (1 - Se):",
          min = 1,
          max = 50,
          value = 15,
          width='100%'
        )
    ),
    div(style = "float: right; width: 48%", # this number can be any percentage you need
        sliderInput(
          "FPR",
          "False Positive Rate: Percentage of true non-patients with test scores above the intersection (1 - Sp):",
          min = 1,
          max = 50,
          value = 15,
          width = '100%'
        )
    ),    
    textOutput("overlap"),
    p(
      "The point of intersection is equal to the optimal dichotomous threshold 
      at which the sum of Sensitivity and Specificity (Se + Sp) is maximized,
      and the sum of errors (FNR + FPR) is minimized. With the sliders 
      the percentages of dichotomous errors and the quality of the test are manipulated."
    ),
    checkboxInput(
      "combineSliders",
      "Combine the FNR slider with the FPR slider ",
      TRUE
    ),

    div(style = "float: left; width: 48%", # this number can be any percentage you need
        sliderInput(
          "negpost.crit",
          "Criterion for negative post-test probability:",
          min = 0.01,
          max = .20,
          value = .10,
          step=.01,
          width='100%'
        )
    ),
    div(style = "float: right; width: 48%", # this number can be any percentage you need
        sliderInput(
          "pospost.crit",
          "Criterion for positive post-test probability:",
          min = .80,
          max = .99,
          value = .90,
          step=.01,
          width = '100%'
        )
    ), 
    radioButtons(
      "outputTable",
      label = h5(
        "Output Table:"
      ),
      choices = list(
        "Sensitivity and Specificity table" = 1,
        "Predictive Values table" = 2
      ),
      inline = T,
      selected = 1
    )
  ),
  
  # Show a plot of the generated distribution
  mainPanel(
    plotOutput("distPlot", height='400px'),
    tableOutput('table'),
    p("* Using the intersection as dichotomous cutpoint."),
    h3(
      "Medical Decision Methods: Grey Zone Determination for Binormal 
       Distributed Test Scores"
    ),
    tabsetPanel(
      tabPanel(
        "Introduction",
        h2('Determination of a Grey zone'),
        p(
          "This Shiny App provides a hands-on demonstration of
          the", span("Grey Zone method for bi-normal distributions", style =
          "color:blue"), ", used for trichotomizing medical test results. In this
          simulation, tests are simulated with a given optimal dichotomous
          False Negatives Rate (1 - Sensitivity) and False Positive Rate (1 -
          Specificity), allowing to compare the trichotomized results with the
          optimal dichotomous results. Tests are simulated with higher
            test scores for indicating the presence of the disease."
          ),
        p(
          "The Grey zone method defines a range of test scores with an
          insufficient post-test probability to distinguishes patients with
          and without the targeted disease (Grey zone), and determines which
          test scores offer a higher or lower post-test probability than the
          pre-set criteria that the targeted disease is present. Only the test
          scores outside the Grey zone are used for a positive or negative
          classification. The patients with test scores inside the Grey zone
          can be retested or further developments can be awaited. These three
          zones offer therefore a trichotomization of the test scores. This
          Shiny App shows advantages and disadvantages of this evaluation
          method of test scores. It is intended to help you to decide for or
          against the use of this method for your medical test."
        ), 
        p(
          "The Grey zone method is developed for epidemiological purposes and
          determines two thresholds, one for distinguishing patients whose
          test score have a higher positive post-test probability than the
          preset value, and one for distinguishing patients whose test score
          has a lower negative post-test probability than the preset value.
          The Grey Zone is the interval of test scores that are inconclusive
          and offer insufficient grounds to distinguish patients with or without
          the targeted disease."
        ),
        h4('Pre-test probability and populations'),
        p(
          'In epidemiological studies, the question is often how you can
          measure disease and health outcomes in a population. Using a test
          that indicates the presence or absence of a disease, a relevant
          question is that of the probability 
          that a group of persons selected from the general population have
          the disease or not, given their test scores. In that case, a suitable pre-test probability is
          the prevalence of the disease in the population.'
          ),
       p(
       'When a clinical population is involved, for instance a group of
          patients who are referred to the Neurology departments of hospitals,
          the prevalence of a disease (for instance Mild Cognitive
          Impairment), may differ strongly from hospital to
          hospital and may differ over time. Clearly, the prevalence of the disease in this
          pre-selected clinical population is higher than in the general
          population. Choosing a suitable pre-test probability may be
          challenging. A pre-test probability of .5 can be suitable,
          reflecting uncertainty about the true prevalence in that clinical
          population.'
          ),
       p(
       'In clinical diagnostic situations, a single test score is obtained and
       the question is whether the patient is selected from the
       sub-population with or from the sub-population without the disease. A
       pre-test probability of .5 is most often suitable, reflecting uncertainty
       about the presence or absence of the disease.'
       ),
       h4('Other Tabs'),
       
        tags$ul(
          tags$li(
            "Hands-on. Look at the Hands-on tab for some examples and their discussion."
          ),
          tags$li(
            "Control panel. Use this tab for explanation of the use of the controls."
          ),
          tags$li("Background. Explanation of the method."),
          tags$li(
            "Conclusion and References. My conclusions (but please draw your own) 
            and References."
          ),
        ),
      ),  # tabPanel
      
      tabPanel(
        "Hands-on",
        h2(
          "Hands-on demonstration of the bi-normal Grey zone method"
        ),
        p(
          # "Set the criterion for the post-test probability to .1 and .9 (start
          # values) and check the check box to combine the two sliders. With the
          # sliders you manipulate the overlap between the two distributions and
          # the strength of the simulated test. A lower overlap directly defines
          # a stronger test. For convenience, the", span("AUC statistic", style
          # = "color:blue"), "is also presented, which is a commonly used
          # estimate of the strength of the test."
        ),
        tags$ul(
          tags$li(
            "1. Move the sliders to values larger than 15%. Observe that the
            Grey zone grows when the test has more overlap and becomes weaker.
            The Grey Interval has a desirable, low post-test probability of .5
            (see Predictive Values Table), which shows that these test scores
            badly distinguish patients from non-patients. Please note that the
            Grey Zone grows larger and larger when the test has too much
            overlap and is therefore of low quality. When the sliders hit 37,
            the size of the Grey Zone is 100% and covers all scores. The Grey
            Zone method is not to be used for weak tests."
          ),
          tags$li(
            '2. Make the test even stronger by reducing the percentages
            further. When moving the FNR slider to values below 10, the two
            borders are switched (black to the right and red to the left) and
            the Grey Zone has a negative size. These stronger tests have a
            better optimal dichotomous result and a Grey zone with a negative
            size'), tags$em('reduces'), ('the result to the criterion. Clearly, such a
            Grey Zone should not be used. Use a higher criterion or be
            satisfied with the dichotomous results.'
          ),
          tags$li(
            "3. Make the test stronger by reducing the overlap between the two
            distributions of the simulated test scores to two times 10%. With
            the criterion of .9 for the positive post probability, the left
            and right border of the Grey Zone coincide with the intersection
            and the Grey Zone does not exist (has a size of 0). The optimal
            dichotomous results are equal. For such a strong test, a higher
            criterion is needed."
          ),
          tags$li(
            "4. Uncheck the check box for combining the two sliders. Keep the
            pre-test probability of .5 and the criteria for the post-test
            probabilities .1 and .9. Now we can generate instances of
            distributions with unequal variance. Set the False Negative Rate
            (FNR) to 15% and the False Positive Rate to 25%. The variance of
            the distribution of true patients is now reduced to .73 and is
            lower than the variance of the non-patients of 1. Please note the
            existence of two thresholds for the positive post-test
            probability. The smallest Grey Zone is used (red vertical line)
            and the secondary threshold (dashed red line) is ignored. The
            existence of two thresholds indicate problematic test scores where
            higher scores are more indicative of absence than presence of the
            disease, lowering the positive post-test probability.", 'Also,
            note that the right border walks away from the intersection and,
            consequently, the Grey Zone includes a part of the distribution
            where the probability of presence of the disease is considerably
            larger than absence of the disease: the post-test probability of
            the Grey Zone is .66 (Predictive Values table). This may be
            considered as not desirable and you may (if possible) want to make
            the Grey Zone smaller. You can do this by selecting a lower value
            for the positive post-test criterion (less variance, less
            restrictive) or a lower value for the negative post-test criterion
            (more variance, more restrictive) to obtain a smaller Grey Zone.
            Change for instance the value for the positive post-test criterion
            to .84. In that case, the post-test predictive value of the Grey
            Zone is .56 (Predictive Values table), which may be considered as
            more desirable, because less patients with the disease are
            excluded, while in this smaller Grey Zone there is only an about
            equal chance to be selected from either distribution.'
            ),
          tags$li(
          "5. Move the FPR slider to 47. Now, the positive post-test
          probability takes a deep dive and does not reach the criterion
          anymore. The desired Grey Zone is not available. Again, you have to
          set a lower criterion here."
          ),
          tags$li(
            "6. Move the FPR slider back to 15 and the post-test probabilities
            on .1 and .9. Move pre-test probability slider slowly to the left.
            This makes the distribution of patients with the disease present
            gradually smaller than the distribution of patients with the
            disease absent. The optimal dichotomous results for Sp and Se
            remains the same: both .85. Also, AUC decreases very slowly.
            Please observe that the Grey Zone walks to the left. Given a test
            with given Sensitivity and Specificity, changing the pre-test
            probability also changes the relative variances of the
            distributions. (Reversely, a test may not have the same Sensitivity
            and Specificity when the relative size of the two
            distributions change. Not demonstrated here.)",
            p("Set the outcomes
            table to Predictive Values. Set the pre-test probability to 39%.
            Please observe that the lower border of the Grey zone is about equal
            to the intersection, while the right border is above the mean of the
            red distribution of patients with the disease present. A relatively
            large part of the patients with the disease fall in the Grey Zone.
            Please note that post-test probability for the lower, negative
            non-Grey interval is .10 and for the upper, positive non-Grey
            interval is .9, exactly as requested by setting the criterion
            values. The Grey zone interval itself has a high post-test
            probability of .73, which may be undesired.")
          ), 
          tags$li(
            "7. Now, set the False Negative Rate to 34%. With a
            pretest-probability of 39 and a FPR of 15, two distributions with
            equal variance are obtained. Now, the Grey Zone is about equally
            distributed around the intersection, with a desirable low
            post-test probability of .49 for the Grey zone."
          ),
        ),
        p(
          "In this simulation, tests are simulated with a given optimal
          dichotomous False Negatives Rate (1 - Sensitivity) and False
          Positive Rate (1 - Specificity), allowing to compare the
          trichotomized results with the optimal dichotomous results. In that
          case, there are two ways to create distributions with unequal
          variance: by choosing different values for FPR and FNR, or by
          setting the pre-probability to a value other than .5, creating a
          relatively smaller or larger distribution of patients with the
          disease. A point of attention is whether the Grey zone is about
          equally distributed around the intersection, which provides the best
          results. When Grey Zone is not around the intersection or does not
          contain the intersection, the post-test predictive value of the test
          scores within the Grey Zone may be undesirable high. One might lower
          or increase the criterion for either the negative post probability
          criterion or the positive post probability criterion, to obtain
          better results where the post-test probability is low for the Grey
          Zone. It is undesirable to not use test scores that have a good
          post-test predictability."
          )
      ), # tabPanel
      
      tabPanel(
        "Control panel",
        h2('Using the dash board'),
        h4('Setting the pre-test probability'),
        p(
          'The grey panel offers a dash board where the user can simulate many
          different tests with bi-normal distributions. The true presence or
          absence of the disease is assumed to be known and is determined with
          superior means, called a "gold standard". The two bi-normal
          distributions form a mixture of patients without (black) and with
          (red) the targeted disease. The pre-test probability is the
          probability of selecting a patient with the disease. When a patient
          is selected from the general population, the pre-test probability
          equals the prevalence of the disease. In a typical clinical
          diagnostic situation, the test is selected to test for a disease
          which presence is probable, given the symptoms of the patient. The
          clinical population is then a population of all patients with the
          given symptoms, who may have or may not have the disease for which
          will be tested. Often, such a clinical population is not well
          described and knowledge of the true presence of the disease, given
          the symptoms is lacking. Therefore, a pre-test probability of .5 may
          be chosen, an equal chance that the patient is chosen from either of
          the two sub-populations. One could say that this reflects the
          pre-test uncertainty whether the patient, given the symptoms, has
          the illness or not.'),
        h4('Setting the False Negative Rate (FNR) and False Positive Rate (FPR'),
        p(
          'The tests differ in quality by their overlap of the scores for the
          two groups. The overlap is chosen with two sliders: The upper slider
          sets the percentage of patients with test scores above the point of
          intersection (blue line in the left graph), indicating the amount of
          False Negatives when the optimal dichotomous cut-point would be
          used. The lower slider sets the percentage of patients which are
          known to not have the disease and who have test scores below the
          point of intersection, indicating the amount of False Negatives for
          the dichotomous optimal cut-point. This allows us to compare the
          results of the trichotomization method with the optimal dichotomous
          method. The FNR equals 1 - Sensitivity and the FPR equals 1 -
          Specificity. The two values specify the weighted overlap between the
          distributions, and the weight is set by the pre-test probability.',
          p(
            "The total overlap between the two distributions defines directly
          the strength of the simulated test, weighted by the pre-test
          probability. The overlap is the percentage of all patients' test
          scores that are part of both the two distributions. The smaller the
          overlap, the stronger the test. In this way, a large amount of
          bi-normal tests of varying strengths can be simulated. The strength
          of the test is directly determined by the overlap of the
          distributions of test scores: a test is stronger when the overlap is
          smaller. "
          ),
        ),
        h4('Combining the FNR slider with the FPR slider'),
        p(
          'A check box allows the combinations of the two sliders. When the two
          sliders are combined, the variance remains equal for the two
          distributions. Unchecking makes it possible to use the two sliders
          separately, allowing the two distributions to have different
          variance.'
        ),
        h4('Criteria for the negative and positive post-test probability'),
        p(
          'In this application to demonstrate the working of the Grey zone
          method, the criteria for the negative and positive post-test
          probability can be set. The post-test probability is the probability
          of the target disorder after a diagnostic test result is known. The
          criterion for the negative post-test probability can be set to
          values between .05 and 0.2, where a lower value indicates a more
          restrictive value. The negative post test probability is the
          probability of the target disorder for patients who receive a
          negative decision based on the test result. The criterion for the
          positive post-test probability can be set to values between .8 and
          0.95, where a higher value indicates a more restrictive value. The
          positive post test probability is the probability of the target
          disorder for patients who receive a positive decision based on the
          test result.'
        ),
        h4('Setting Output table'),
        p('The output can be set to a Sensitivity and Specificity table, or to
        a table with Predictive Values (including the post-test
        probability).')

      ), # tabPanel
      tabPanel("Output Tables",
       h2('Sensitivity and Specificity table'),
       p(
         'The first row provides the Se and Sp for the dichotomized results.
         The results are optimized by using the intersection as dichotomous
         threshold, which maximizes the sum of Se + Sp and minimizes the sum
         of FNR + FPR. The dichotomization uses all scores. The calculations
         of Se and Sp are as usual and are not dependent on the pre-test probability.'
       ),
       p(
         'The second row provides the Sensitivity and Specificity solely for
         the test scores within the Grey Zone, with the intersection as
         cut-point. This can only be calculated when the point of intersection
         falls within the Grey Zone. Desirable values for both Se and Sp of
         the Grey zone are close to .5, indicating that the test is neither
         Sensitive nor Specific when using these scores. These Se and Sp are
         dependent on the pre-test probability, as the  Grey zone may contain
         a larger part from one distribution and a smaller part of the other.
         Also, the intersection may divide the Grey Zone in unequal parts.'
       ),
       p(
         'The third row provides the Sensitivity and Specificity solely for
         the test scores outside the Grey Zone. This can only be calculated
         when the point of intersection falls within the Grey Zone. Desirable
         values for both Se and Sp of the Grey zone are values higher than the
         dichotomous results, indicating that the test is more Sensitive and/or
         Specific when using these scores. These Se and Sp are dependent on
         the pre-test probability, as the  Grey zone may contain a larger part
         from one distribution and a smaller part of the other.'
       ),
       h2('Predictive Values table'),
       p(
         'Predictive values can be calculated for every interval of test
         scores. The post-test predictive probabilities are always equal to
         the Positive Predictive Values (PPV). The Standardized Positive
         Predictive Value is equal to the PPV, only when the pre-test
         probability is .5. The calculated results for the post-test
         predictive probability of the lower negative non-grey interval of
         test scores is always equal to the set criterion for the negative
         post-test probability. Similarly is the equality between the
         post-test predictive probability of the upper positive non-grey
         interval of test scores and the set criterion for the positive
         post-test probability. The post-test probabilities of the lower
         and upper non-Grey Zone are always equal to the pre-set criteria for
         the negative and positive post-test probability. Desirable values are
         better than the dichotomous positive post-test probability (and 1 -
         this value for the set negative criterion value).'
         ),
       p(
         'The first row shows the predictive values for the dichotomized
         results in the usual way with the intersection as threshold. Here,
         the post predictive value is the positive post-test predictive value;
         the negative post-test predictive value is 1 - the positive post-test
         predictive value (not shown). Please note that the values of the
         Standardized Negative Predictive Value (SNPV) and of the Standardized
         Positive Predictive Value (SPPV) are independent of the Pre-test
         probability. All other values are dependent on the Pre-test
         probability.'
       ),
      p(
      'The second row shows the predictive values for solely the test-scores
      within the Grey Zone. Please note the SNPV and SPPV are not independent
      on the pre-test probability as the position of the Grey Zone may concern
      a relative larger part of one of the distributions.'
      ),
      p(
      'The third row shows the predictive values for solely the test-scores of
      the negative non-grey interval below the Grey Zone. Please note the SNPV
      and SPPV are not independent on the pre-test probability. '
      ),
      p(
      'The fourth row shows the predictive values for solely the test-scores
      of the positive non-grey interval above the Grey Zone. Please note the
      SNPV and SPPV are not independent on the pre-test probability. '
      ),
      ), # tabPanel
      
      tabPanel(
        "Background",
        h2("Background Information"),
        p(
            'This method is proposed by Coste et al. (2003, 2006). This method
            uses all possible test scores as dichotomous thresholds to
            calculate Se, Sp, positive and negative likelihood ratios and
            post-test probabilities. The post-test probabilities are calculated for
            the accumulated densities of the test scores and indicate the
            levels of seriousness of the disease for test scores above and
            below all possible dichotomous thresholds.'),
        p(
           'The criterion is a required degree of closeness of post-test
           probabilities to 1 or 0. The default criterion values are .05 and
           .95 for respectively a negative and positive classification, which
           may be sufficient for use by clinicians or Public Health
           professionals for a first classification whether a target disease
           may be present or not (Coste et al., 2003). Assuming that a test
           distinguishes patients from non-patients, the pre-test probability
           sets the lower limit of the possible positive post-test
           probabilities and the upper limit of the possible negative post-tets
           probabilities. The upper limit of the possible positive post-test
           probabilities and the lower limit of the negative post-test
           probabilities are dependent on the test quality. Therefore, possible
           criteria for the post-test probabilities are limited by the pre-test
           probability and the quality of the test.'
          ),
        p(
          'A test, intended for confirming the presence or absence of a
          disease, can be evaluated using two groups: a group of true patients
          who have the disease and a group of patients who truly do not have
          the disease (shortly called non-patients). The plot shows the two
          simulated densities of the the test scores of the two groups. Higher
          test scores indicate the presence of the disease. The ', span( 'Grey
          zone method', style = "color:blue"), 'determines the borders of a
          Grey zone for these test scores. It is assumed here that the test
          scores are normally distributed. This Grey Zone indicates a range of
          test scores that have sufficient post test probabilities to
          distinguish the two groups of patients. '
        ),
        
        h2('Test score Densities'),
        p(
          "The distributions shown in the plot show the densities of the
          obtained simulated test scores, which are assumed to be normal
          distributed. The test scores of the non-patients are standard normal
          distributed, with mean of 0 and standard deviation of 1: N(0, 1).
          The distribution of the true patients can vary widely. For a given
          test score, the difference of the two densities indicates from which
          of the two groups of patients the test score is most likely drawn.
          When the sliders are combined, both simulated distributions have a
          standard deviation of 1 and only the means differs. The application
          starts with a distribution of test scores with a mean of 2.07 and a
          standard deviation of 1 (N(2.07, 1))."
        ),
        
        h2('Using Densities'),
        p(
          'In the Grey Zone method, post-test probabilities are used for the
          determination whether a testcore (or a small range of test scores)
          is sampled from the group of patients or from the group of
          non-patients. The Grey zone method is intended for classification of
          patients, of who it is suspected that they may have the targeted
          disease. At the start, the pretest probability is set to .50. In
          that case, before testing  it is considered equally likely that
          the patient may or may not have the disease.'
        ),
        
        h2("Trichotomization versus dichotomization"),
        p(
          "The classic techniques for evaluating tests for medical decision
          make use of dichotomization of the test scores. Implicitly, all
          test scores are considered as equally valid, and all test scores are
          used for a positive or negative classification. The R package
          UncertainInterval offers several trichotomization methods, that 
          identify test scores that are insufficiently valid and are better
          not used for classification. The idea behind these trichotomization
          methods is that patients are better served when these invalid test
          scores are not used for a positive or negative classification. These
          patients are better off with retesting, if possible with a better
          test, or awaiting further developments. In some medical fields the
          uncertainty of a diagnostic outcome is followed up by techniques
          such as watchful waiting and active surveillance. These more
          cautious approaches are especially relevant when possible treatments
          can have serious side effects."
        ),
      ), # tabPanel
      
      tabPanel(
        "Conclusion and References",
        h2('In conclusion'),
        p(
          'The largest advantage of using trichotomization methods is
          the identification of patients with uncertain, inconclusive test
          scores. This allows for following a more cautious path for these
          patients.'),
        p(
        'The simulations show that the Grey Zone method for bi-normal data may
        offer a better positive or negative classification of patients, with a
        sufficient post-test probability, when solely the test scores outside
        the Grey Zone are used. To reach this goal, criteria for the post-test
        probability have to be chosen, relative to the quality of the test.
        The method works best when the variances of the distributions of both
        the sub-populations of patients are about equal. When the variance is
        unequal, a part of scores may be placed in the Grey Zone, that in fact
        offer a reasonable distinction between the two sub-populations and therefore a
        relatively high pre-test probability of presence of the disease in
        the Grey Zone. In that case, the criterion for the sub-population with
        the higher variance has to be lowered.'
        ),
        p(
          'Assuming that a test distinguishes patients from non-patients, the
          pre-test probability sets the lower limit of the possible positive
          post-test probabilities and the upper limit of the possible negative
          post-tets probabilities. The upper limit of the possible positive
          post-test probabilities and the lower limit of the negative
          post-test probabilities are dependent on the test quality. Therefore,
          possible criteria for the post-test probabilities are limited by the
          pre-test probability and the quality of the test.'
          ),
        p(
          'The main disadvantage of this method is that different criteria
          have to be applied, dependent on the the pre-test probability and
          the quality of the test. Important for the qualities of the test are
          the variances of the test scores in the two sub-populations of
          patients with and without the targeted disease. Setting suitable
          criteria is therefore especially a problem when clinical populations
          are involved for which the distributions are not exactly known.'
          ),
        h2('References'),
        p(
          'Coste, J., Jourdain, P., & Pouchot, J. (2006). A gray zone assigned
          to inconclusive results of quantitative diagnostic tests:
          application to the use of brain natriuretic peptide for diagnosis of
          heart failure in acute dyspneic patients. Clinical Chemistry,
          52(12), 2229-2235.'
        ),
        p(
          'Coste, J., & Pouchot, J. (2003). A grey zone for quantitative
          diagnostic and screening tests. International Journal of
          Epidemiology, 32(2), 304-313.'
        ),
        p(
          'Landsheer, J. A. (2020). Impact of the Prevalence of Cognitive
             Impairment on the Accuracy of the Montreal Cognitive Assessment:
             The Advantage of Using two MoCA Thresholds to Identify Error-prone
             Test Scores. Alzheimer Disease & Associated Disorders.',
          tags$a(href = "https://dx.doi.org/10.1097/WAD.0000000000000365",
                 "https://dx.doi.org/10.1097/WAD.0000000000000365")
        ),
        p(
          'Landsheer J. A. (2020). CRAN - Package UncertainInterval 0.6.0.',
          tags$a(
            href = " https://CRAN.R-project.org/package=UncertainInterval",
            " https://CRAN.R-project.org/package=UncertainInterval"
          )
        ),
        p(
          'Youden, W. J. (1950). Index for rating diagnostic tests.
           Cancer, 3(1), 32–35.',
          tags$a(
            href = "https://doi.org/10.1002/1097-0142(1950)3:1<32::AID-CNCR2820030106>3.0.CO;2-3",
            "https://doi.org/10.1002/1097-0142(1950)3:1<32::AID-CNCR2820030106>3.0.CO;2-3"
          )
        ),
        
      ) # tabPanel
  ) # tabsetpanel
  ) # mainpanel
  ) # sidebarLayout
) # fluidpage
