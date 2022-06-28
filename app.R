#
# Simulator of P-Value Multiple Hypothesis Correction
# (SIMPLYCORRECT)
# This is a Shiny web application for visualizing three p-value
# correction methods--Bonferroni, Benjamini-Hochberg, and
# permutation FDR control--and examining the costs and benefits
# of each.

# 
# Run runBeforePublishing.R before publishing the app.
#

library(shiny)
library(shinyjs)
source("SIMPLYCORRECT_functions.R")

ui <- navbarPage(
  "SIMPLYCORRECT",
  
  #### Welcome Tab ####
  tabPanel(
    "Welcome",
    
    # Application title
    titlePanel("Simulator of P-Value Multiple Hypothesis Correction (SIMPLYCORRECT)"),
    
    fluidRow(column(width = 12, h3("Welcome!"))),
    
    fluidRow(column(width = 12, h4("To get started quickly, click a tab at the top."))),
    
    fluidRow(column(width = 12, downloadLink("downloadData2", "CLICK HERE TO DOWNLOAD MANUAL!"))),
    
    br(),
    
    fluidRow(column(width = 12, p("This simulator is a tool for visualizing the costs and effects of correcting P-values in quantitative omics experiments."))),
    
    fluidRow(column(width = 12, p("As an educational tool, SIMPLYCORRECT can illustrate the importance of P-value correction in various theoretical circumstances. It can also illustrate how to increase the sensitivity of quantitative omics experiments by increasing sample sizes and decreasing variability."))),
    
    fluidRow(column(width = 12, p("As a research tool, if estimates of experiment parameters are known, e.g., from a pilot experiment, these parameter estimates can be plugged into the model to estimate the sensitivity and specificity of full-scale studies on the same system as the pilot. However, use with extreme caution, as these simple models assume analyte levels are mutually independent, which is usually not true. Dependence is highly common in quantitative proteomics, transcriptomics, and metabolomics and will affect results."))),
    
    fluidRow(column(width = 12, downloadLink("downloadData", "For a detailed instruction manual with explanations of parameters, click here."))),
    
    br()
  ),
  
  #### Tab 1: No changes ####
  tabPanel(
    "1: No Analytes Change",
    
    # Application title
    titlePanel("Simulator of P-Value Multiple Hypothesis Correction (SIMPLYCORRECT)"),
    
    fluidRow(column(width = 12, h3("    Model 1: No Analytes Change"))),
    
    fluidRow(column(width = 12, downloadLink("downloadData3", "CLICK HERE TO DOWNLOAD MANUAL!"))),
    
    fluidRow(column(width = 12, checkboxInput("toggleTables", "Show result tables", value = T))),
    
    br(),
    
    # Side bar page layout
    sidebarLayout(
      
      # Side bar
      sidebarPanel(
        width = 3,
        
        #Style
        tags$head(
          tags$style(HTML('#goButton{background-color:#88BBFF;
                          border: 2px solid black;
                          width: 100%;
                          font-weight: bold}'),
                     HTML('#BHResults{margin: auto;
                          padding: 30px}'),
                     HTML('#unadjResults{margin: auto;
                          padding: 30px}'),
                     HTML('#BonfResults{margin: auto;
                          padding: 30px}'),
                     HTML('#PermResults{margin: auto;
                          padding: 30px}')
          )
        ),
        
        # Choose which adjustments to perform
        strong("Correction methods: "),
        checkboxInput("correctBH", "Benjamini-Hochberg", value = T),
        checkboxInput("correctBonf", "Bonferroni", value = F),
        checkboxInput("correctPerm", "Permutation FDR", value = F),
        
        # Standard options
        numericInput("m", "Number of analytes:", value = 400),
        numericInput("nc", "Number of samples in control group:", value = 3),
        numericInput("nt", "Number of samples in test group:", value = 3),
        numericInput("threshold", "Significance threshold:", value = 0.05, step = 0.01),
        numericInput("expVar", "Experimental variability within analyte (C.V.):", value = 0.01, step = 0.01),
        
        # Enable advanced options checkbox
        checkboxInput("advanced", "Show advanced settings", value = F),
        
        # Advanced options
        conditionalPanel(
          condition = "input.advanced == true",
          
          # Number of permutations if doing perm FDR
          conditionalPanel(
            condition = "input.correctPerm == true",
            numericInput("nperms", "Number of permutations:", value = 100)
          ),
          
          # Avg quantity and SD, if box is checked
          numericInput("mu", "Average true quantity across analytes:", value = 22),
          numericInput("sigma", "Standard deviation across analytes:", value = 3),
          # Numeric input for setting the seed
          numericInput("mySeed", "Set seed:", value =1)
        ),
        
        # Make plot button
        actionButton("goButton", "Go!")
        
      ),
      
      # Main panel
      # Show a volcano plot with unadjusted P-values
      # and another for each P-value adjustment method
      mainPanel(
        
        # Top row of plots
        fluidRow(
          
          # Unadj plot
          column(6,plotOutput(outputId = "unadjVolcano")),
          
          # BH plot
          conditionalPanel(
            condition = "input.correctBH == true",
            column(6,plotOutput(outputId = "BHVolcano"))
          )
        ),
        
        conditionalPanel(
          condition = "input.toggleTables == true",
          
          # Unadj table
          fluidRow(
            column(6,tableOutput(outputId = "unadjResults")),
            
            # BH table
            conditionalPanel(
              condition = "input.correctBH == true",
              column(6,tableOutput(outputId = "BHResults"))
            )
          )
          
        ), br(),
        
        # Second row of plots
        fluidRow(
          
          # Bonferroni plot
          conditionalPanel(
            condition = "input.correctBonf == true",
            column(6,plotOutput(outputId = "BonfVolcano"))
          ),
          
          # Perm FDR plot
          conditionalPanel(
            condition = "input.correctPerm == true",
            column(6,plotOutput(outputId = "PermVolcano"))
          )
        ),
        
        conditionalPanel(
          condition = "input.toggleTables == true",
          
          fluidRow(
            # Bonferroni table
            conditionalPanel(
              condition = "input.correctBonf == true",
              column(6,tableOutput(outputId = "BonfResults"))
            ),
            
            # Perm FDR table
            conditionalPanel(
              condition = "input.correctPerm == true",
              column(6,tableOutput(outputId = "PermResults"))
            )
          )
          
        )
      )
    ),
    br(),
    br()
  ),
  
  #### Tab 2: All analytes change ####
  tabPanel(
    "2: All Analytes Change",
    
    # Application title
    titlePanel("Simulator of P-Value Multiple Hypothesis Correction (SIMPLYCORRECT)"),
    
    fluidRow(column(width = 12, h3("    Model 2: All Analytes Change"))),
    
    fluidRow(column(width = 12, downloadLink("downloadData4", "CLICK HERE TO DOWNLOAD MANUAL!"))),
    
    fluidRow(column(width = 12, checkboxInput("toggleTables2", "Show result tables", value = T))),
    
    br(),
    
    # Side bar page layout
    sidebarLayout(
      
      # Side bar
      sidebarPanel(
        width = 3,
        
        #Style
        tags$head(
          tags$style(HTML('#goButton2{background-color:#88BBFF;
                          border: 2px solid black;
                          width: 100%;
                          font-weight: bold}'),
                     HTML('#BHResults2{margin: auto;
                          padding: 30px}'),
                     HTML('#unadjResults2{margin: auto;
                          padding: 30px}'),
                     HTML('#BonfResults2{margin: auto;
                          padding: 30px}'),
                     HTML('#PermResults2{margin: auto;
                          padding: 30px}')
          ),
          tags$style(
            type="text/css", "#effect2Label{color: blue}"
          )
        ),
        
        # Choose which adjustments to perform
        strong("Correction methods: "),
        checkboxInput("correctBH2", "Benjamini-Hochberg", value = T),
        checkboxInput("correctBonf2", "Bonferroni", value = F),
        checkboxInput("correctPerm2", "Permutation FDR", value = F),
        
        # Tab 2 option
        div( id = "effect2Label",
             numericInput("effect2", "True effect:", value = 0.5)
        ),
        
        # Standard options
        numericInput("m2", "Number of analytes:", value = 400),
        numericInput("nc2", "Number of samples in control group:", value = 4),
        numericInput("nt2", "Number of samples in test group:", value = 4),
        numericInput("threshold2", "Significance threshold:", value = 0.05, step = 0.01),
        numericInput("expVar2", "Experimental variability within analyte (C.V.):", value = 0.01, step = 0.01),
        
        # Enable advanced options checkbox
        checkboxInput("advanced2", "Show advanced settings", value = F),
        
        # Advanced options
        conditionalPanel(
          condition = "input.advanced2 == true",
          
          # Number of permutations if doing perm FDR
          conditionalPanel(
            condition = "input.correctPerm2 == true",
            numericInput("nperms2", "Number of permutations:", value = 100)
          ),
          
          # Avg quantity and SD, if box is checked
          numericInput("mu2", "Average true quantity across analytes:", value = 22),
          numericInput("sigma2", "Standard deviation across analytes:", value = 3),
          # Numeric input for setting the seed
          numericInput("mySeed2", "Set seed:", value =1)
        ),
        
        # Make plot button
        actionButton("goButton2", "Go!")
        
      ),
      
      # Main panel
      # Show a volcano plot with unadjusted P-values
      # and another for each P-value adjustment method
      mainPanel(
        
        # Top row of plots
        fluidRow(
          
          # Unadj plot
          column(6,plotOutput(outputId = "unadjVolcano2")),
          
          # BH plot
          conditionalPanel(
            condition = "input.correctBH2 == true",
            column(6,plotOutput(outputId = "BHVolcano2"))
          )
        ),
        
        conditionalPanel(
          condition = "input.toggleTables2 == true",
          
          # Unadj table
          fluidRow(
            column(6,tableOutput(outputId = "unadjResults2")),
            
            # BH table
            conditionalPanel(
              condition = "input.correctBH2 == true",
              column(6,tableOutput(outputId = "BHResults2"))
            )
          )
          
        ), br(),
        
        # Second row
        fluidRow(
          
          # Bonferroni plot
          conditionalPanel(
            condition = "input.correctBonf2 == true",
            column(6,plotOutput(outputId = "BonfVolcano2"))
          ),
          
          # Perm FDR plot
          conditionalPanel(
            condition = "input.correctPerm2 == true",
            column(6,plotOutput(outputId = "PermVolcano2"))
          )
        ),
        
        conditionalPanel(
          condition = "input.toggleTables2 == true",
          
          fluidRow(
            # Bonferroni table
            conditionalPanel(
              condition = "input.correctBonf2 == true",
              column(6,tableOutput(outputId = "BonfResults2"))
            ),
            
            # Perm FDR table
            conditionalPanel(
              condition = "input.correctPerm2 == true",
              column(6,tableOutput(outputId = "PermResults2"))
            )
          )
          
        )
      )
    ),
    
    br(),
    br()
    
  ),
  
  #### Tab 3: Realistic model ####
  tabPanel(
    "3: Realistic Model",
    
    # Application title
    titlePanel("Simulator of P-Value Multiple Hypothesis Correction (SIMPLYCORRECT)"),
    
    fluidRow(column(width = 12, h3("    Model 3: Realistic Model"))),
    
    fluidRow(column(width = 12, downloadLink("downloadData5", "CLICK HERE TO DOWNLOAD MANUAL!"))),
    
    fluidRow(column(width = 12, checkboxInput("toggleTables3", "Show result tables", value = F))),
    
    br(),
    
    # Side bar page layout
    sidebarLayout(
      
      # Side bar
      sidebarPanel(
        width = 3,
        
        #Style
        tags$head(
          tags$style(HTML('#goButton3{background-color:#88BBFF;
                          border: 2px solid black;
                          width: 100%;
                          font-weight: bold}'),
                     HTML('#BHResults3{margin: auto;
                          padding: 30px}'),
                     HTML('#unadjResults3{margin: auto;
                          padding: 30px}'),
                     HTML('#BonfResults3{margin: auto;
                          padding: 30px}'),
                     HTML('#PermResults3{margin: auto;
                          padding: 30px}')
          ),
          tags$style(
            type="text/css", "#effect3Label{color: blue}",
            "#bioVar3Label{color: blue}",
            "#expVar3Label{color: blue}"
          )
        ),
        
        # Choose which adjustments to perform
        strong("Correction methods: "),
        checkboxInput("correctBH3", "Benjamini-Hochberg", value = T),
        checkboxInput("correctBonf3", "Bonferroni", value = T),
        checkboxInput("correctPerm3", "Permutation FDR", value = T),
        
        # Tab 3 option
        div( id = "effect3Label",
             numericInput("effect3", "Typical effect:", value = 1)
        ),
        div( id = "bioVar3Label",
             numericInput("bioVar3", "Typical biological variability (C.V.):", value = 0.05)
        ),
        div(id = "expVar3Label",
            numericInput("expVar3", "Technical variability (C.V.):", value = 0.03, step = 0.01)
        ),
        
        # Standard options
        numericInput("m3", "Number of analytes:", value = 800),
        numericInput("nc3", "Number of samples in control group:", value = 10),
        numericInput("nt3", "Number of samples in test group:", value = 10),
        numericInput("threshold3", "Significance threshold:", value = 0.05, step = 0.01),
        
        # Enable advanced options checkbox
        checkboxInput("advanced3", "Show advanced settings", value = F),
        
        # Advanced options
        conditionalPanel(
          condition = "input.advanced3 == true",
          
          # Number of permutations if doing perm FDR
          conditionalPanel(
            condition = "input.correctPerm3 == true",
            numericInput("nperms3", "Number of permutations:", value = 100)
          ),
          
          # Avg quantity and SD, if box is checked
          numericInput("mu3", "Average true quantity across analytes:", value = 22),
          numericInput("sigma3", "Standard deviation across analytes:", value = 3),
          # Numeric input for setting the seed
          numericInput("mySeed3", "Set seed:", value =1),
          
          # Biological variability distribution parameters
          hr(),
          em(h4("Biological variability distribution parameters")),
          numericInput("bioVarK3", "Integer shape parameter (k):", value = 4, step = 1),
          numericInput("bioVarTheta3", "Scale parameter (theta):", value = 0.5, step = 0.1)
        ),
        
        # Make plot button
        actionButton("goButton3", "Go!")
        
      ),
      
      # Main panel
      # Show a volcano plot with unadjusted P-values
      # and another for each P-value adjustment method
      mainPanel(
        
        # Top row of plots
        fluidRow(
          
          # Unadj plot
          column(6,plotOutput(outputId = "unadjVolcano3")),
          
          # BH plot
          conditionalPanel(
            condition = "input.correctBH3 == true",
            column(6,plotOutput(outputId = "BHVolcano3"))
          )
        ),
        
        conditionalPanel(
          condition = "input.toggleTables3 == true",
          
          # Unadj table
          fluidRow(
            column(6,tableOutput(outputId = "unadjResults3")),
            
            # BH table
            conditionalPanel(
              condition = "input.correctBH3 == true",
              column(6,tableOutput(outputId = "BHResults3"))
            )
          )
          
        ), br(),
        
        # Second row
        fluidRow(
          
          # Bonferroni plot
          conditionalPanel(
            condition = "input.correctBonf3 == true",
            column(6,plotOutput(outputId = "BonfVolcano3"))
          ),
          
          # Perm FDR plot
          conditionalPanel(
            condition = "input.correctPerm3 == true",
            column(6,plotOutput(outputId = "PermVolcano3"))
          )
        ),
        
        conditionalPanel(
          condition = "input.toggleTables3 == true",
          
          fluidRow(
            # Bonferroni table
            conditionalPanel(
              condition = "input.correctBonf3 == true",
              column(6,tableOutput(outputId = "BonfResults3"))
            ),
            
            # Perm FDR table
            conditionalPanel(
              condition = "input.correctPerm3 == true",
              column(6,tableOutput(outputId = "PermResults3"))
            )
          )
          
        )
      )
    ),
    
    br(),
    br()
    
  )
  
  
)
# Define server logic required to draw a histogram
server <- function(input, output) {
  
  #### Welcome tab ####
  output$downloadData <- downloadHandler(
    filename = "22-03-25_SRS_SIMPLYCORRECT-Manual.pdf",
    content = function(file) {
      file.copy("www/22-03-25_SRS_SIMPLYCORRECT-Manual.pdf", file)
    }
  )
  
  output$downloadData2 <- downloadHandler(
    filename = "22-03-25_SRS_SIMPLYCORRECT-Manual.pdf",
    content = function(file) {
      file.copy("www/22-03-25_SRS_SIMPLYCORRECT-Manual.pdf", file)
    }
  )
  
  output$downloadData3 <- downloadHandler(
    filename = "22-03-25_SRS_SIMPLYCORRECT-Manual.pdf",
    content = function(file) {
      file.copy("www/22-03-25_SRS_SIMPLYCORRECT-Manual.pdf", file)
    }
  )
  
  output$downloadData4 <- downloadHandler(
    filename = "22-03-25_SRS_SIMPLYCORRECT-Manual.pdf",
    content = function(file) {
      file.copy("www/22-03-25_SRS_SIMPLYCORRECT-Manual.pdf", file)
    }
  )
  
  output$downloadData5 <- downloadHandler(
    filename = "22-03-25_SRS_SIMPLYCORRECT-Manual.pdf",
    content = function(file) {
      file.copy("www/22-03-25_SRS_SIMPLYCORRECT-Manual.pdf", file)
    }
  )
  
  #### Tab 1 functions ####  
  makeUnadjData = eventReactive(input$goButton, ignoreInit = F, ignoreNULL = F, {
    simulateData(seed = input$mySeed,
                 m = input$m,
                 avg = input$mu,
                 sd = input$sigma,
                 nc = input$nc,
                 nt = input$nt,
                 expVar = input$expVar,
                 threshold = input$threshold)
  })
  
  output$unadjVolcano = renderPlot({
    dataAndColors = makeUnadjData()
    unadjData = dataAndColors[[1]]
    myColors = dataAndColors[[2]]
    threshold = dataAndColors[[4]]
    volcanoPlot(unadjData, myColors, threshold)
  })
  
  output$unadjResults = renderTable({
    resultTable = makeUnadjData()
    resultTable[[3]]
  })
  
  makeBHData = eventReactive(input$goButton, ignoreInit = F, ignoreNULL = F, {
    simulateData(seed = input$mySeed,
                 m = input$m,
                 avg = input$mu,
                 sd = input$sigma,
                 nc = input$nc,
                 nt = input$nt,
                 expVar = input$expVar,
                 threshold = input$threshold,
                 adjMeth = "BH")
  })
  
  output$BHVolcano = renderPlot({
    BHList = makeBHData()
    BHData = BHList[[1]]
    myColors = BHList[[2]]
    BHThreshold = BHList[[4]]
    volcanoPlot(BHData, myColors, BHThreshold, "Benjamini-Hochberg Correction")
  })
  
  output$BHResults = renderTable({
    resultTable = makeBHData()
    resultTable[[3]]
  })
  
  makeBonfData = eventReactive(input$goButton, ignoreInit = F, ignoreNULL = T, {
    simulateData(seed = input$mySeed,
                 m = input$m,
                 avg = input$mu,
                 sd = input$sigma,
                 nc = input$nc,
                 nt = input$nt,
                 expVar = input$expVar,
                 threshold = input$threshold,
                 adjMeth = "Bonf")
  })
  
  output$BonfVolcano = renderPlot({
    BonfList = makeBonfData()
    BonfData = BonfList[[1]]
    myColors = BonfList[[2]]
    BonfThreshold = BonfList[[4]]
    volcanoPlot(BonfData, myColors, BonfThreshold, "Bonferroni Correction")
  })
  
  output$BonfResults = renderTable({
    resultTable = makeBonfData()
    resultTable[[3]]
  })
  
  makePermData = eventReactive(input$goButton, ignoreInit = F, ignoreNULL = T, {
    simulateData(seed = input$mySeed,
                 m = input$m,
                 avg = input$mu,
                 sd = input$sigma,
                 nc = input$nc,
                 nt = input$nt,
                 expVar = input$expVar,
                 threshold = input$threshold,
                 adjMeth = "Perm")
  })
  
  output$PermVolcano = renderPlot({
    PermList = makePermData()
    PermData = PermList[[1]]
    myColors = PermList[[2]]
    PermThreshold = PermList[[4]]
    volcanoPlot(PermData, myColors, PermThreshold, "Permutation FDR Correction")
  })
  
  output$PermResults = renderTable({
    resultTable = makePermData()
    resultTable[[3]]
  })
  
  #### Tab 2 functions ####
  
  makeUnadjData2 = eventReactive(input$goButton2, ignoreInit = F, ignoreNULL = F, {
    if (input$effect2 == 0) {
      
      simulateData(seed = input$mySeed2,
                   m = input$m2,
                   avg = input$mu2,
                   sd = input$sigma2,
                   nc = input$nc2,
                   nt = input$nt2,
                   expVar = input$expVar2,
                   threshold = input$threshold2,
                   model = "NoChange",
                   effectMagnitude = input$effect2)
      
    } else {
      
      simulateData(seed = input$mySeed2,
                   m = input$m2,
                   avg = input$mu2,
                   sd = input$sigma2,
                   nc = input$nc2,
                   nt = input$nt2,
                   expVar = input$expVar2,
                   threshold = input$threshold2,
                   model = "AllChange",
                   effectMagnitude = input$effect2)
    }
  })
  
  output$unadjVolcano2 = renderPlot({
    dataAndColors = makeUnadjData2()
    unadjData = dataAndColors[[1]]
    myColors = dataAndColors[[2]]
    threshold = dataAndColors[[4]]
    volcanoPlot(unadjData, myColors, threshold)
  })
  
  output$unadjResults2 = renderTable({
    resultTable = makeUnadjData2()
    resultTable[[3]]
  })
  
  makeBHData2 = eventReactive(input$goButton2, ignoreInit = F, ignoreNULL = F, {
    
    if (input$effect2 == 0) {
      
      simulateData(seed = input$mySeed2,
                   m = input$m2,
                   avg = input$mu2,
                   sd = input$sigma2,
                   nc = input$nc2,
                   nt = input$nt2,
                   expVar = input$expVar2,
                   threshold = input$threshold2,
                   adjMeth = "BH",
                   model = "NoChange",
                   effectMagnitude = input$effect2)
      
    } else {  
      
      simulateData(seed = input$mySeed2,
                   m = input$m2,
                   avg = input$mu2,
                   sd = input$sigma2,
                   nc = input$nc2,
                   nt = input$nt2,
                   expVar = input$expVar2,
                   threshold = input$threshold2,
                   adjMeth = "BH",
                   model = "AllChange",
                   effectMagnitude = input$effect2)
    }
  })
  
  output$BHVolcano2 = renderPlot({
    BHList = makeBHData2()
    BHData = BHList[[1]]
    myColors = BHList[[2]]
    BHThreshold = BHList[[4]]
    volcanoPlot(BHData, myColors, BHThreshold, "Benjamini-Hochberg Correction")
  })
  
  output$BHResults2 = renderTable({
    resultTable = makeBHData2()
    resultTable[[3]]
  })
  
  makeBonfData2 = eventReactive(input$goButton2, ignoreInit = F, ignoreNULL = T, {
    if (input$effect2 == 0) {
      
      simulateData(seed = input$mySeed2,
                   m = input$m2,
                   avg = input$mu2,
                   sd = input$sigma2,
                   nc = input$nc2,
                   nt = input$nt2,
                   expVar = input$expVar2,
                   threshold = input$threshold2,
                   adjMeth = "Bonf",
                   model = "NoChange",
                   effectMagnitude = input$effect2)
      
    } else {
      
      simulateData(seed = input$mySeed2,
                   m = input$m2,
                   avg = input$mu2,
                   sd = input$sigma2,
                   nc = input$nc2,
                   nt = input$nt2,
                   expVar = input$expVar2,
                   threshold = input$threshold2,
                   adjMeth = "Bonf",
                   model = "AllChange",
                   effectMagnitude = input$effect2)
    }
  })
  
  output$BonfVolcano2 = renderPlot({
    BonfList = makeBonfData2()
    BonfData = BonfList[[1]]
    myColors = BonfList[[2]]
    BonfThreshold = BonfList[[4]]
    volcanoPlot(BonfData, myColors, BonfThreshold, "Bonferroni Correction")
  })
  
  output$BonfResults2 = renderTable({
    resultTable = makeBonfData2()
    resultTable[[3]]
  })
  
  makePermData2 = eventReactive(input$goButton2, ignoreInit = F, ignoreNULL = T, {
    if (input$effect2 == 0) {
      
      simulateData(seed = input$mySeed2,
                   m = input$m2,
                   avg = input$mu2,
                   sd = input$sigma2,
                   nc = input$nc2,
                   nt = input$nt2,
                   expVar = input$expVar2,
                   threshold = input$threshold2,
                   adjMeth = "Perm",
                   model = "NoChange",
                   effectMagnitude = input$effect2)
      
    } else {
      
      simulateData(seed = input$mySeed2,
                   m = input$m2,
                   avg = input$mu2,
                   sd = input$sigma2,
                   nc = input$nc2,
                   nt = input$nt2,
                   expVar = input$expVar2,
                   threshold = input$threshold2,
                   adjMeth = "Perm",
                   model = "AllChange",
                   effectMagnitude = input$effect2)
    }
  })
  
  output$PermVolcano2 = renderPlot({
    PermList = makePermData2()
    PermData = PermList[[1]]
    myColors = PermList[[2]]
    PermThreshold = PermList[[4]]
    volcanoPlot(PermData, myColors, PermThreshold, "Permutation FDR Correction")
  })
  
  output$PermResults2 = renderTable({
    resultTable = makePermData2()
    resultTable[[3]]
  })
  
  
  
  #### Tab 3 functions ####
  
  makeUnadjData3 = eventReactive(input$goButton3, ignoreInit = F, ignoreNULL = F, {
    simulateData(seed = input$mySeed3,
                 m = input$m3,
                 avg = input$mu3,
                 sd = input$sigma3,
                 nc = input$nc3,
                 nt = input$nt3,
                 expVar = input$expVar3,
                 threshold = input$threshold3,
                 model = "Realistic",
                 effectMagnitude = input$effect3,
                 effectSD = input$effect3,
                 bioVar = input$bioVar3,
                 bioVGamma = c(input$bioVarK3, input$bioVarTheta3)
    )
  })
  
  output$unadjVolcano3 = renderPlot({
    dataAndColors = makeUnadjData3()
    unadjData = dataAndColors[[1]]
    myColors = dataAndColors[[2]]
    threshold = dataAndColors[[4]]
    volcanoPlot(unadjData, myColors, threshold)
  })
  
  output$unadjResults3 = renderTable({
    resultTable = makeUnadjData3()
    resultTable[[3]]
  })
  
  makeBHData3 = eventReactive(input$goButton3, ignoreInit = F, ignoreNULL = F, {
    simulateData(seed = input$mySeed3,
                 m = input$m3,
                 avg = input$mu3,
                 sd = input$sigma3,
                 nc = input$nc3,
                 nt = input$nt3,
                 expVar = input$expVar3,
                 threshold = input$threshold3,
                 adjMeth = "BH",
                 model = "Realistic",
                 effectMagnitude = input$effect3,
                 effectSD = input$effect3,
                 bioVar = input$bioVar3,
                 bioVGamma = c(input$bioVarK3, input$bioVarTheta3))
  })
  
  output$BHVolcano3 = renderPlot({
    BHList = makeBHData3()
    BHData = BHList[[1]]
    myColors = BHList[[2]]
    BHThreshold = BHList[[4]]
    volcanoPlot(BHData, myColors, BHThreshold, "Benjamini-Hochberg Correction")
  })
  
  output$BHResults3 = renderTable({
    resultTable = makeBHData3()
    resultTable[[3]]
  })
  
  makeBonfData3 = eventReactive(input$goButton3, ignoreInit = F, ignoreNULL = F, {
    simulateData(seed = input$mySeed3,
                 m = input$m3,
                 avg = input$mu3,
                 sd = input$sigma3,
                 nc = input$nc3,
                 nt = input$nt3,
                 expVar = input$expVar3,
                 threshold = input$threshold3,
                 adjMeth = "Bonf",
                 model = "Realistic",
                 effectMagnitude = input$effect3,
                 effectSD = input$effect3,
                 bioVar = input$bioVar3,
                 bioVGamma = c(input$bioVarK3, input$bioVarTheta3))
  })
  
  output$BonfVolcano3 = renderPlot({
    BonfList = makeBonfData3()
    BonfData = BonfList[[1]]
    myColors = BonfList[[2]]
    BonfThreshold = BonfList[[4]]
    volcanoPlot(BonfData, myColors, BonfThreshold, "Bonferroni Correction")
  })
  
  output$BonfResults3 = renderTable({
    resultTable = makeBonfData3()
    resultTable[[3]]
  })
  
  makePermData3 = eventReactive(input$goButton3, ignoreInit = F, ignoreNULL = F, {
    simulateData(seed = input$mySeed3,
                 m = input$m3,
                 avg = input$mu3,
                 sd = input$sigma3,
                 nc = input$nc3,
                 nt = input$nt3,
                 expVar = input$expVar3,
                 threshold = input$threshold3,
                 adjMeth = "Perm",
                 model = "Realistic",
                 effectMagnitude = input$effect3,
                 effectSD = input$effect3,
                 bioVar = input$bioVar3,
                 bioVGamma = c(input$bioVarK3, input$bioVarTheta3))
  })
  
  output$PermVolcano3 = renderPlot({
    PermList = makePermData3()
    PermData = PermList[[1]]
    myColors = PermList[[2]]
    PermThreshold = PermList[[4]]
    volcanoPlot(PermData, myColors, PermThreshold, "Permutation FDR Correction")
  })
  
  output$PermResults3 = renderTable({
    resultTable = makePermData3()
    resultTable[[3]]
  })
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)

