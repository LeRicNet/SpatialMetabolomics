#' SpatialMetabolics Interactive Explorer
#'
#' Shiny application for interactive exploration of SpatialMetabolic objects
#'
#' @author Your Name

library(shiny)
library(shinydashboard)
library(SpatialMetabolics)
library(ggplot2)
library(plotly)
library(DT)

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "SpatialMetabolics Explorer"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("home")),
      menuItem("Spatial Plots", tabName = "spatial", icon = icon("map")),
      menuItem("Pathway Analysis", tabName = "pathways", icon = icon("project-diagram")),
      menuItem("Differential Expression", tabName = "de", icon = icon("chart-bar")),
      menuItem("Quality Control", tabName = "qc", icon = icon("check-circle"))
    ),

    hr(),

    # File input
    fileInput("file", "Load SpatialMetabolic Object",
              accept = c(".rds", ".RDS")),

    # Or use example data
    actionButton("useExample", "Use Example Data",
                 icon = icon("database"),
                 class = "btn-primary btn-block")
  ),

  dashboardBody(
    tabItems(
      # Overview tab
      tabItem(tabName = "overview",
              h2("Dataset Overview"),
              fluidRow(
                valueBoxOutput("nSamples"),
                valueBoxOutput("nSpots"),
                valueBoxOutput("nGenes")
              ),
              fluidRow(
                box(
                  title = "Sample Information",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 6,
                  DT::dataTableOutput("sampleTable")
                ),
                box(
                  title = "Pathway Coverage",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 6,
                  plotOutput("pathwayCoverage")
                )
              )
      ),

      # Spatial plots tab
      tabItem(tabName = "spatial",
              h2("Spatial Visualization"),
              fluidRow(
                box(
                  title = "Controls",
                  status = "primary",
                  width = 3,
                  selectInput("spatialFeature", "Feature to plot:",
                              choices = NULL),
                  radioButtons("featureType", "Feature type:",
                               choices = c("Pathway Score" = "score",
                                           "Gene Expression" = "gene")),
                  selectInput("spatialSample", "Sample:",
                              choices = NULL,
                              multiple = TRUE),
                  selectInput("colorScale", "Color scale:",
                              choices = c("viridis", "plasma", "inferno", "magma"),
                              selected = "viridis"),
                  sliderInput("pointSize", "Point size:",
                              min = 0.1, max = 3, value = 1, step = 0.1),
                  actionButton("updateSpatial", "Update Plot",
                               icon = icon("refresh"),
                               class = "btn-primary")
                ),
                box(
                  title = "Spatial Plot",
                  status = "primary",
                  width = 9,
                  plotlyOutput("spatialPlot", height = "600px")
                )
              )
      ),

      # Pathway analysis tab
      tabItem(tabName = "pathways",
              h2("Pathway Analysis"),
              fluidRow(
                box(
                  title = "Pathway Scores",
                  status = "primary",
                  width = 6,
                  plotlyOutput("pathwayBoxplot")
                ),
                box(
                  title = "Pathway Correlations",
                  status = "primary",
                  width = 6,
                  plotlyOutput("pathwayHeatmap")
                )
              ),
              fluidRow(
                box(
                  title = "Metabolic State",
                  status = "primary",
                  width = 12,
                  column(4,
                         selectInput("pathway1", "Pathway 1:",
                                     choices = NULL)
                  ),
                  column(4,
                         selectInput("pathway2", "Pathway 2:",
                                     choices = NULL)
                  ),
                  column(4,
                         selectInput("scatterGroup", "Color by:",
                                     choices = NULL)
                  ),
                  plotlyOutput("metabolicScatter", height = "400px")
                )
              )
      ),

      # Differential expression tab
      tabItem(tabName = "de",
              h2("Differential Expression"),
              fluidRow(
                box(
                  title = "Analysis Settings",
                  status = "primary",
                  width = 3,
                  selectInput("deCondition", "Condition column:",
                              choices = NULL),
                  selectInput("refGroup", "Reference group:",
                              choices = NULL),
                  selectInput("testGroup", "Test group:",
                              choices = NULL),
                  radioButtons("deFeatures", "Features to test:",
                               choices = c("Pathway Scores" = "scores",
                                           "All Genes" = "genes")),
                  actionButton("runDE", "Run Analysis",
                               icon = icon("play"),
                               class = "btn-success")
                ),
                box(
                  title = "Results",
                  status = "primary",
                  width = 9,
                  conditionalPanel(
                    condition = "output.deComplete",
                    tabsetPanel(
                      tabPanel("Volcano Plot",
                               plotlyOutput("volcanoPlot", height = "500px")
                      ),
                      tabPanel("Results Table",
                               DT::dataTableOutput("deTable")
                      ),
                      tabPanel("Top Features",
                               plotOutput("topFeatures", height = "600px")
                      )
                    )
                  )
                )
              )
      ),

      # QC tab
      tabItem(tabName = "qc",
              h2("Quality Control"),
              fluidRow(
                box(
                  title = "QC Metrics",
                  status = "primary",
                  width = 12,
                  plotlyOutput("qcPlots", height = "600px")
                )
              ),
              fluidRow(
                box(
                  title = "Sample Correlation",
                  status = "primary",
                  width = 6,
                  plotlyOutput("sampleCorrelation")
                ),
                box(
                  title = "PCA",
                  status = "primary",
                  width = 6,
                  plotlyOutput("pcaPlot")
                )
              )
      )
    )
  )
)

# Define server
server <- function(input, output, session) {

  # Reactive values
  values <- reactiveValues(
    spm = NULL,
    de_results = NULL
  )

  # Load data
  observeEvent(input$file, {
    req(input$file)

    withProgress(message = "Loading data...", {
      values$spm <- loadSpatialMetabolic(input$file$datapath)
    })

    updateUI(values$spm)
  })

  # Use example data
  observeEvent(input$useExample, {
    withProgress(message = "Loading example data...", {
      # Create example data
      values$spm <- create_test_data(n_genes = 500, n_spots = 200)
      values$spm <- normalizeSpatial(values$spm, verbose = FALSE)
      values$spm <- calculateQCMetrics(values$spm)

      # Calculate some pathways
      pathways <- getAllMetabolicPathways("mouse")[1:3]
      values$spm <- calculateMetabolicScores(values$spm, pathways, verbose = FALSE)
    })

    updateUI(values$spm)
  })

  # Update UI elements when data is loaded
  updateUI <- function(spm) {
    # Update feature selections
    pathways <- rownames(metabolicScores(spm))
    genes <- head(rownames(spm), 100)  # Limit to first 100 genes

    updateSelectInput(session, "spatialFeature",
                      choices = c(pathways, genes))

    updateSelectInput(session, "pathway1", choices = pathways)
    updateSelectInput(session, "pathway2", choices = pathways)

    # Update sample selections
    samples <- unique(colData(spm)$sample_id)
    updateSelectInput(session, "spatialSample",
                      choices = samples,
                      selected = samples[1])

    # Update condition selections
    meta_cols <- colnames(colData(spm))
    factor_cols <- meta_cols[sapply(colData(spm), is.factor)]

    updateSelectInput(session, "deCondition", choices = factor_cols)
    updateSelectInput(session, "scatterGroup", choices = factor_cols)
  }

  # Update group selections
  observeEvent(input$deCondition, {
    req(values$spm, input$deCondition)

    groups <- unique(colData(values$spm)[[input$deCondition]])
    updateSelectInput(session, "refGroup", choices = groups)
    updateSelectInput(session, "testGroup", choices = groups)
  })

  # Overview outputs
  output$nSamples <- renderValueBox({
    valueBox(
      value = ifelse(is.null(values$spm), 0,
                     length(unique(colData(values$spm)$sample_id))),
      subtitle = "Samples",
      icon = icon("vials"),
      color = "blue"
    )
  })

  output$nSpots <- renderValueBox({
    valueBox(
      value = ifelse(is.null(values$spm), 0, ncol(values$spm)),
      subtitle = "Spots",
      icon = icon("circle"),
      color = "green"
    )
  })

  output$nGenes <- renderValueBox({
    valueBox(
      value = ifelse(is.null(values$spm), 0, nrow(values$spm)),
      subtitle = "Genes",
      icon = icon("dna"),
      color = "yellow"
    )
  })

  # Sample table
  output$sampleTable <- DT::renderDataTable({
    req(values$spm)

    sample_info <- as.data.frame(colData(values$spm)) %>%
      group_by(sample_id) %>%
      summarise(
        n_spots = n(),
        mean_counts = mean(nCount_RNA),
        mean_features = mean(nFeature_RNA),
        .groups = "drop"
      )

    DT::datatable(sample_info, options = list(pageLength = 10))
  })

  # Pathway coverage plot
  output$pathwayCoverage <- renderPlot({
    req(values$spm)

    if (nrow(metabolicScores(values$spm)) > 0) {
      pathway_list <- metabolicPathways(values$spm)

      coverage_data <- data.frame(
        pathway = names(pathway_list),
        total = sapply(pathway_list, length),
        present = sapply(pathway_list, function(x) sum(x %in% rownames(values$spm)))
      )
      coverage_data$coverage <- coverage_data$present / coverage_data$total * 100

      ggplot(coverage_data, aes(x = pathway, y = coverage)) +
        geom_col(fill = "steelblue") +
        geom_text(aes(label = paste0(round(coverage, 1), "%")),
                  vjust = -0.5) +
        theme_minimal() +
        labs(x = "", y = "Coverage (%)")
    }
  })

  # Spatial plot
  output$spatialPlot <- renderPlotly({
    req(input$updateSpatial, values$spm)

    isolate({
      p <- plotSpatialMetabolicScore(
        values$spm,
        pathway = input$spatialFeature,
        samples = input$spatialSample,
        color_scale = input$colorScale,
        point_size = input$pointSize
      )

      ggplotly(p)
    })
  })

  # Pathway boxplot
  output$pathwayBoxplot <- renderPlotly({
    req(values$spm)

    if (nrow(metabolicScores(values$spm)) > 0) {
      p <- plotMetabolicSummary(values$spm, plot_type = "box")
      ggplotly(p)
    }
  })

  # Run differential expression
  observeEvent(input$runDE, {
    req(values$spm, input$deCondition, input$refGroup, input$testGroup)

    withProgress(message = "Running differential expression...", {
      if (input$deFeatures == "scores") {
        features <- NULL
      } else {
        features == rownames(values$spm)
      }

      values$de_results <- compareMetabolicStates(
        values$spm,
        condition = input$deCondition,
        ref_group = input$refGroup,
        test_group = input$testGroup,
        features = features,
        verbose = FALSE
      )
    })
  })

  # DE complete flag
  output$deComplete <- reactive({
    !is.null(values$de_results)
  })
  outputOptions(output, "deComplete", suspendWhenHidden = FALSE)

  # Volcano plot
  output$volcanoPlot <- renderPlotly({
    req(values$de_results)

    p <- plotVolcano(values$de_results, label_top = 10)
    ggplotly(p)
  })

  # DE table
  output$deTable <- DT::renderDataTable({
    req(values$de_results)

    DT::datatable(
      values$de_results,
      options = list(
        pageLength = 20,
        scrollX = TRUE
      )
    ) %>%
      formatRound(columns = c("log2FC", "pvalue", "adj_pvalue"), digits = 3)
  })
}

# Run app
shinyApp(ui = ui, server = server)
