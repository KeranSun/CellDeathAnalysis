#' =============================================================================
#' Shiny Interactive Application for CellDeathAnalysis
#' =============================================================================

#' Launch CellDeathAnalysis Shiny App
#'
#' Launch an interactive Shiny application for cell death pathway analysis.
#'
#' @param max_upload_size Maximum file upload size in MB. Default is 100.
#' @param ... Additional arguments passed to shiny::runApp.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch the app
#' launch_death_app()
#'
#' # Launch with custom settings
#' launch_death_app(max_upload_size = 200, port = 3838)
#' }
#'
launch_death_app <- function(max_upload_size = 100, ...) {
 
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Install with: install.packages('shiny')")
  }
 
  if (!requireNamespace("shinydashboard", quietly = TRUE)) {
    stop("Package 'shinydashboard' is required. Install with: install.packages('shinydashboard')")
  }
 
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("Package 'DT' is required. Install with: install.packages('DT')")
  }
 
  # Set upload size
  options(shiny.maxRequestSize = max_upload_size * 1024^2)
 
  # Create app
  app <- shiny::shinyApp(
    ui = .death_app_ui(),
    server = .death_app_server
  )
 
  shiny::runApp(app, ...)
}


#' Internal: Shiny UI
#' @keywords internal
.death_app_ui <- function() {
 
  shinydashboard::dashboardPage(
    skin = "blue",
   
    # Header
    shinydashboard::dashboardHeader(
      title = "CellDeathAnalysis",
      titleWidth = 250
    ),
   
    # Sidebar
    shinydashboard::dashboardSidebar(
      width = 250,
      shinydashboard::sidebarMenu(
        id = "tabs",
        shinydashboard::menuItem("Home", tabName = "home",
                                  icon = shiny::icon("home")),
        shinydashboard::menuItem("Data Upload", tabName = "upload",
                                  icon = shiny::icon("upload")),
        shinydashboard::menuItem("Gene Sets", tabName = "genesets",
                                  icon = shiny::icon("dna")),
        shinydashboard::menuItem("Pathway Scoring", tabName = "scoring",
                                  icon = shiny::icon("calculator")),
        shinydashboard::menuItem("Visualization", tabName = "visualization",
                                  icon = shiny::icon("chart-bar")),
        shinydashboard::menuItem("Survival Analysis", tabName = "survival",
                                  icon = shiny::icon("heartbeat")),
        shinydashboard::menuItem("Enrichment", tabName = "enrichment",
                                  icon = shiny::icon("search")),
        shinydashboard::menuItem("Machine Learning", tabName = "ml",
                                  icon = shiny::icon("brain")),
        shinydashboard::menuItem("Export Results", tabName = "export",
                                  icon = shiny::icon("download")),
        shinydashboard::menuItem("Help", tabName = "help",
                                  icon = shiny::icon("question-circle"))
      )
    ),
   
    # Body
    shinydashboard::dashboardBody(
      # CSS
      shiny::tags$head(
        shiny::tags$style(shiny::HTML("
          .content-wrapper { background-color: #f4f6f9; }
          .box { box-shadow: 0 1px 3px rgba(0,0,0,0.12); }
          .info-box { min-height: 90px; }
          .nav-tabs-custom > .tab-content { padding: 15px; }
        "))
      ),
     
      shinydashboard::tabItems(
       
        # Home Tab
        shinydashboard::tabItem(
          tabName = "home",
          shiny::fluidRow(
            shiny::column(12,
              shinydashboard::box(
                title = "Welcome to CellDeathAnalysis",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                shiny::h4("Comprehensive Analysis of Cell Death Pathways"),
                shiny::p("This interactive application allows you to analyze cell death
                          pathway-related gene expression in your transcriptomic data."),
                shiny::hr(),
                shiny::h5("Supported Cell Death Types:"),
                shiny::p("Ferroptosis, Cuproptosis, Disulfidptosis, Pyroptosis, Necroptosis,
                          Apoptosis, Autophagy, PANoptosis, NETosis, Parthanatos, Entosis,
                          Oxeiptosis, Alkaliptosis, LDCD"),
                shiny::hr(),
                shiny::h5("Quick Start:"),
                shiny::tags$ol(
                  shiny::tags$li("Upload your expression data (CSV/TSV/RDS)"),
                  shiny::tags$li("Calculate pathway scores"),
                  shiny::tags$li("Visualize and analyze results"),
                  shiny::tags$li("Export results")
                )
              )
            )
          ),
          shiny::fluidRow(
            shinydashboard::infoBox(
              "Pathways", 14, icon = shiny::icon("dna"),
              color = "blue", width = 3
            ),
            shinydashboard::infoBox(
              "Genes", "500+", icon = shiny::icon("database"),
              color = "green", width = 3
            ),
            shinydashboard::infoBox(
              "Methods", 6, icon = shiny::icon("calculator"),
              color = "yellow", width = 3
            ),
            shinydashboard::infoBox(
              "Version", "0.3.0", icon = shiny::icon("code-branch"),
              color = "red", width = 3
            )
          )
        ),
       
        # Upload Tab
        shinydashboard::tabItem(
          tabName = "upload",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Upload Expression Data",
              status = "primary",
              solidHeader = TRUE,
              width = 6,
              shiny::fileInput("expr_file", "Expression Matrix (CSV/TSV/RDS)",
                               accept = c(".csv", ".tsv", ".txt", ".rds")),
              shiny::helpText("Rows = genes, Columns = samples"),
              shiny::checkboxInput("expr_header", "First row is header", TRUE),
              shiny::checkboxInput("expr_rownames", "First column is row names", TRUE),
              shiny::hr(),
              shiny::actionButton("load_example", "Load Example Data",
                                  icon = shiny::icon("database"),
                                  class = "btn-info")
            ),
            shinydashboard::box(
              title = "Upload Clinical Data (Optional)",
              status = "info",
              solidHeader = TRUE,
              width = 6,
              shiny::fileInput("clinical_file", "Clinical Data (CSV/TSV)",
                               accept = c(".csv", ".tsv", ".txt")),
              shiny::selectInput("sample_col", "Sample ID Column", choices = NULL),
              shiny::selectInput("group_col", "Group Column", choices = NULL),
              shiny::selectInput("time_col", "Survival Time Column", choices = NULL),
              shiny::selectInput("status_col", "Survival Status Column", choices = NULL)
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Data Summary",
              status = "success",
              solidHeader = TRUE,
              width = 12,
              shiny::verbatimTextOutput("data_summary")
            )
          )
        ),
       
        # Gene Sets Tab
        shinydashboard::tabItem(
          tabName = "genesets",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Cell Death Pathways",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              DT::DTOutput("pathway_table")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Pathway Details",
              status = "info",
              solidHeader = TRUE,
              width = 6,
              shiny::selectInput("selected_pathway", "Select Pathway",
                                 choices = NULL),
              shiny::verbatimTextOutput("pathway_info")
            ),
            shinydashboard::box(
              title = "Gene List",
              status = "info",
              solidHeader = TRUE,
              width = 6,
              DT::DTOutput("gene_list")
            )
          )
        ),
       
        # Scoring Tab
        shinydashboard::tabItem(
          tabName = "scoring",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Scoring Parameters",
              status = "primary",
              solidHeader = TRUE,
              width = 4,
              shiny::selectInput("score_method", "Scoring Method",
                                 choices = c("Z-score" = "zscore",
                                             "Mean" = "mean",
                                             "Median" = "median",
                                             "ssGSEA" = "ssgsea",
                                             "GSVA" = "gsva",
                                             "AUCell" = "aucell")),
              shiny::checkboxGroupInput("score_pathways", "Select Pathways",
                                        choices = NULL),
              shiny::actionButton("btn_select_all", "Select All"),
              shiny::actionButton("btn_deselect_all", "Deselect All"),
              shiny::hr(),
              shiny::actionButton("run_scoring", "Calculate Scores",
                                  icon = shiny::icon("play"),
                                  class = "btn-success btn-lg")
            ),
            shinydashboard::box(
              title = "Score Results",
              status = "success",
              solidHeader = TRUE,
              width = 8,
              shiny::verbatimTextOutput("score_summary"),
              DT::DTOutput("score_table")
            )
          )
        ),
       
        # Visualization Tab
        shinydashboard::tabItem(
          tabName = "visualization",
          shiny::fluidRow(
            shinydashboard::tabBox(
              title = "Visualization",
              width = 12,
              shiny::tabPanel(
                "Heatmap",
                shiny::plotOutput("plot_heatmap", height = "500px"),
                shiny::downloadButton("dl_heatmap", "Download")
              ),
              shiny::tabPanel(
                "Boxplot",
                shiny::selectInput("boxplot_pathways", "Select Pathways",
                                   choices = NULL, multiple = TRUE),
                shiny::plotOutput("plot_boxplot", height = "500px"),
                shiny::downloadButton("dl_boxplot", "Download")
              ),
              shiny::tabPanel(
                "Radar Chart",
                shiny::plotOutput("plot_radar", height = "500px"),
                shiny::downloadButton("dl_radar", "Download")
              ),
              shiny::tabPanel(
                "Correlation",
                shiny::plotOutput("plot_correlation", height = "500px"),
                shiny::downloadButton("dl_correlation", "Download")
              ),
              shiny::tabPanel(
                "Distribution",
                shiny::selectInput("dist_pathway", "Select Pathway", choices = NULL),
                shiny::selectInput("dist_type", "Plot Type",
                                   choices = c("Density" = "density",
                                               "Histogram" = "histogram",
                                               "Violin" = "violin")),
                shiny::plotOutput("plot_distribution", height = "400px"),
                shiny::downloadButton("dl_distribution", "Download")
              )
            )
          )
        ),
       
        # Survival Tab
        shinydashboard::tabItem(
          tabName = "survival",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Survival Analysis Settings",
              status = "primary",
              solidHeader = TRUE,
              width = 4,
              shiny::selectInput("surv_pathway", "Select Pathway", choices = NULL),
              shiny::selectInput("surv_method", "Grouping Method",
                                 choices = c("Median" = "median",
                                             "Mean" = "mean",
                                             "Quantile" = "quantile")),
              shiny::conditionalPanel(
                condition = "input.surv_method == 'quantile'",
                shiny::sliderInput("surv_quantile", "Quantile", 0.5, 0.9, 0.75, 0.05)
              ),
              shiny::actionButton("run_survival", "Run Analysis",
                                  icon = shiny::icon("play"),
                                  class = "btn-success")
            ),
            shinydashboard::box(
              title = "Kaplan-Meier Curve",
              status = "success",
              solidHeader = TRUE,
              width = 8,
              shiny::plotOutput("plot_km", height = "450px"),
              shiny::downloadButton("dl_km", "Download")
            )
          ),
          shiny::fluidRow(
            shinydashboard::box(
              title = "Survival Statistics",
              status = "info",
              solidHeader = TRUE,
              width = 6,
              shiny::verbatimTextOutput("surv_stats")
            ),
            shinydashboard::box(
              title = "Batch Survival Analysis",
              status = "warning",
              solidHeader = TRUE,
              width = 6,
              shiny::actionButton("run_batch_surv", "Run Batch Analysis",
                                  class = "btn-warning"),
              DT::DTOutput("batch_surv_table")
            )
          )
        ),
       
        # Enrichment Tab
        shinydashboard::tabItem(
          tabName = "enrichment",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Gene Input",
              status = "primary",
              solidHeader = TRUE,
              width = 4,
              shiny::textAreaInput("gene_input", "Enter Gene List (one per line)",
                                   rows = 10,
                                   placeholder = "GPX4\nSLC7A11\nACSL4\n..."),
              shiny::actionButton("run_ora", "Run ORA",
                                  icon = shiny::icon("search"),
                                  class = "btn-success")
            ),
            shinydashboard::box(
              title = "Enrichment Results",
              status = "success",
              solidHeader = TRUE,
              width = 8,
              DT::DTOutput("ora_table"),
              shiny::plotOutput("plot_ora", height = "400px"),
              shiny::downloadButton("dl_ora", "Download Results")
            )
          )
        ),
       
        # ML Tab
        shinydashboard::tabItem(
          tabName = "ml",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Model Settings",
              status = "primary",
              solidHeader = TRUE,
              width = 4,
              shiny::selectInput("ml_outcome", "Outcome Variable", choices = NULL),
              shiny::selectInput("ml_method", "Model Type",
                                 choices = c("Random Forest" = "rf",
                                             "LASSO" = "lasso",
                                             "XGBoost" = "xgboost",
                                             "Ensemble" = "ensemble")),
              shiny::sliderInput("ml_train_ratio", "Training Ratio",
                                 0.5, 0.9, 0.7, 0.05),
              shiny::actionButton("run_ml", "Train Model",
                                  icon = shiny::icon("cogs"),
                                  class = "btn-success")
            ),
            shinydashboard::box(
              title = "Model Performance",
              status = "success",
              solidHeader = TRUE,
              width = 8,
              shiny::verbatimTextOutput("ml_results"),
              shiny::plotOutput("plot_importance", height = "350px")
            )
          )
        ),
       
        # Export Tab
        shinydashboard::tabItem(
          tabName = "export",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Export Results",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              shiny::h4("Download Options"),
              shiny::downloadButton("dl_scores_csv", "Download Scores (CSV)",
                                    class = "btn-info"),
              shiny::downloadButton("dl_scores_rds", "Download Scores (RDS)",
                                    class = "btn-info"),
              shiny::downloadButton("dl_report", "Download Report (HTML)",
                                    class = "btn-success")
            )
          )
        ),
       
        # Help Tab
        shinydashboard::tabItem(
          tabName = "help",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Documentation",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              shiny::h4("Quick Guide"),
              shiny::tags$ol(
                shiny::tags$li(shiny::strong("Data Upload:"),
                               "Upload your expression matrix (genes x samples).
                                Supported formats: CSV, TSV, RDS."),
                shiny::tags$li(shiny::strong("Scoring:"),
                               "Select scoring method and pathways. Z-score is fastest;
                                ssGSEA/GSVA recommended for publication."),
                shiny::tags$li(shiny::strong("Visualization:"),
                               "Create heatmaps, boxplots, radar charts, etc."),
                shiny::tags$li(shiny::strong("Survival:"),
                               "Requires clinical data with time and status columns."),
                shiny::tags$li(shiny::strong("Export:"),
                               "Download results as CSV or RDS files.")
              ),
              shiny::hr(),
              shiny::h4("Contact"),
              shiny::p("Author: Keran Sun"),
              shiny::p("Email: s1214844197@163.com"),
              shiny::p("GitHub: https://github.com/keransun/CellDeathAnalysis")
            )
          )
        )
      )
    )
  )
}


#' Internal: Shiny Server
#' @keywords internal
.death_app_server <- function(input, output, session) {
 
  # Reactive values
  rv <- shiny::reactiveValues(
    expr = NULL,
    clinical = NULL,
    scores = NULL,
    surv_result = NULL,
    ml_model = NULL
  )
 
  # Initialize pathway choices
  shiny::observe({
    pathways <- names(get_death_geneset("all", type = "simple"))
    shiny::updateCheckboxGroupInput(session, "score_pathways",
                                     choices = pathways, selected = pathways)
    shiny::updateSelectInput(session, "selected_pathway", choices = pathways)
    shiny::updateSelectInput(session, "surv_pathway", choices = pathways)
    shiny::updateSelectInput(session, "boxplot_pathways", choices = pathways,
                             selected = pathways[1:4])
    shiny::updateSelectInput(session, "dist_pathway", choices = pathways)
  })
 
  # Pathway table
  output$pathway_table <- DT::renderDT({
    info <- list_death_pathways(detailed = TRUE)
    DT::datatable(info, options = list(pageLength = 15))
  })
 
  # Pathway info
  output$pathway_info <- shiny::renderPrint({
    shiny::req(input$selected_pathway)
    get_pathway_info(input$selected_pathway)
  })
 
  # Gene list
  output$gene_list <- DT::renderDT({
    shiny::req(input$selected_pathway)
    genes <- get_death_geneset(input$selected_pathway, type = "all")
    df <- data.frame(Gene = genes)
    DT::datatable(df, options = list(pageLength = 20))
  })
 
  # Load example data
  shiny::observeEvent(input$load_example, {
    data("example_expr", package = "CellDeathAnalysis", envir = environment())
    data("example_clinical", package = "CellDeathAnalysis", envir = environment())
    rv$expr <- example_expr
    rv$clinical <- example_clinical
   
    shiny::showNotification("Example data loaded!", type = "message")
   
    # Update clinical columns
    cols <- colnames(rv$clinical)
    shiny::updateSelectInput(session, "sample_col", choices = cols, selected = "sample_id")
    shiny::updateSelectInput(session, "group_col", choices = c("None", cols), selected = "group")
    shiny::updateSelectInput(session, "time_col", choices = c("None", cols), selected = "OS_time")
    shiny::updateSelectInput(session, "status_col", choices = c("None", cols), selected = "OS_status")
    shiny::updateSelectInput(session, "ml_outcome", choices = cols)
  })
 
  # Load expression data
  shiny::observeEvent(input$expr_file, {
    file <- input$expr_file
    ext <- tools::file_ext(file$datapath)
   
    if (ext == "rds") {
      rv$expr <- readRDS(file$datapath)
    } else {
      rv$expr <- read.csv(file$datapath,
                           header = input$expr_header,
                           row.names = if(input$expr_rownames) 1 else NULL)
      rv$expr <- as.matrix(rv$expr)
    }
   
    shiny::showNotification("Expression data loaded!", type = "message")
  })
 
  # Load clinical data
  shiny::observeEvent(input$clinical_file, {
    file <- input$clinical_file
    rv$clinical <- read.csv(file$datapath, header = TRUE)
   
    cols <- colnames(rv$clinical)
    shiny::updateSelectInput(session, "sample_col", choices = cols)
    shiny::updateSelectInput(session, "group_col", choices = c("None", cols))
    shiny::updateSelectInput(session, "time_col", choices = c("None", cols))
    shiny::updateSelectInput(session, "status_col", choices = c("None", cols))
    shiny::updateSelectInput(session, "ml_outcome", choices = cols)
   
    shiny::showNotification("Clinical data loaded!", type = "message")
  })
 
  # Data summary
  output$data_summary <- shiny::renderPrint({
    if (is.null(rv$expr)) {
      cat("No expression data loaded.\n")
      cat("Please upload data or click 'Load Example Data'.\n")
    } else {
      cat("Expression Data:\n")
      cat("  Genes:", nrow(rv$expr), "\n")
      cat("  Samples:", ncol(rv$expr), "\n")
      cat("  Range:", round(min(rv$expr), 2), "-", round(max(rv$expr), 2), "\n\n")
     
      if (!is.null(rv$clinical)) {
        cat("Clinical Data:\n")
        cat("  Samples:", nrow(rv$clinical), "\n")
        cat("  Variables:", ncol(rv$clinical), "\n")
      }
    }
  })
 
  # Select all pathways
  shiny::observeEvent(input$btn_select_all, {
    pathways <- names(get_death_geneset("all", type = "simple"))
    shiny::updateCheckboxGroupInput(session, "score_pathways",
                                     selected = pathways)
  })
 
  # Deselect all pathways
  shiny::observeEvent(input$btn_deselect_all, {
    shiny::updateCheckboxGroupInput(session, "score_pathways", selected = character(0))
  })
 
  # Calculate scores
  shiny::observeEvent(input$run_scoring, {
    shiny::req(rv$expr, input$score_pathways)
   
    shiny::withProgress(message = "Calculating scores...", {
      rv$scores <- calculate_death_score(
        rv$expr,
        pathways = input$score_pathways,
        method = input$score_method,
        verbose = FALSE
      )
    })
   
    shiny::showNotification("Scoring complete!", type = "message")
  })
 
  # Score summary
  output$score_summary <- shiny::renderPrint({
    shiny::req(rv$scores)
    cat("Score Matrix:", nrow(rv$scores), "samples x", ncol(rv$scores), "pathways\n\n")
    print(summary(rv$scores[, 1:min(4, ncol(rv$scores))]))
  })
 
  # Score table
  output$score_table <- DT::renderDT({
    shiny::req(rv$scores)
    df <- round(rv$scores, 4)
    df <- cbind(Sample = rownames(df), df)
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
  })
 
  # Get group variable
  get_group <- shiny::reactive({
    if (is.null(rv$clinical) || input$group_col == "None") {
      return(NULL)
    }
    rv$clinical[[input$group_col]]
  })
 
  # Heatmap
  output$plot_heatmap <- shiny::renderPlot({
    shiny::req(rv$scores)
    tryCatch({
      plot_death_heatmap(rv$scores, show_annotation = !is.null(get_group()),
                         annotation = get_group())
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "Heatmap requires ComplexHeatmap or pheatmap package")
    })
  })
 
  # Boxplot
  output$plot_boxplot <- shiny::renderPlot({
    shiny::req(rv$scores, input$boxplot_pathways)
    plot_death_boxplot(rv$scores, group = get_group(),
                       pathways = input$boxplot_pathways)
  })
 
  # Radar
  output$plot_radar <- shiny::renderPlot({
    shiny::req(rv$scores)
    plot_death_radar(rv$scores, group = get_group())
  })
 
  # Correlation
  output$plot_correlation <- shiny::renderPlot({
    shiny::req(rv$scores)
    tryCatch({
      plot_pathway_correlation(rv$scores)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "Correlation plot requires additional packages")
    })
  })
 
  # Distribution
  output$plot_distribution <- shiny::renderPlot({
    shiny::req(rv$scores, input$dist_pathway)
    plot_score_distribution(rv$scores, pathways = input$dist_pathway,
                            type = input$dist_type, group = get_group())
  })
 
  # Run survival analysis
  shiny::observeEvent(input$run_survival, {
    shiny::req(rv$scores, rv$clinical, input$surv_pathway)
    shiny::req(input$time_col != "None", input$status_col != "None")
   
    time <- rv$clinical[[input$time_col]]
    status <- rv$clinical[[input$status_col]]
   
    quantile_val <- if (input$surv_method == "quantile") input$surv_quantile else 0.5
   
    tryCatch({
      rv$surv_result <- death_survival(
        rv$scores,
        time = time,
        status = status,
        pathway = input$surv_pathway,
        method = input$surv_method,
        quantile = quantile_val
      )
      shiny::showNotification("Survival analysis complete!", type = "message")
    }, error = function(e) {
      shiny::showNotification(paste("Error:", e$message), type = "error")
    })
  })
 
  # KM plot
  output$plot_km <- shiny::renderPlot({
    shiny::req(rv$surv_result)
    tryCatch({
      p <- plot_survival(rv$surv_result)
      print(p)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "Survival plot requires survminer package")
    })
  })
 
  # Survival stats
  output$surv_stats <- shiny::renderPrint({
    shiny::req(rv$surv_result)
    print(rv$surv_result)
  })
 
  # ORA analysis
  shiny::observeEvent(input$run_ora, {
    shiny::req(input$gene_input)
   
    genes <- strsplit(input$gene_input, "\n")[[1]]
    genes <- trimws(genes)
    genes <- genes[genes != ""]
   
    if (length(genes) < 2) {
      shiny::showNotification("Please enter at least 2 genes", type = "warning")
      return()
    }
   
    rv$ora_result <- death_enrich_ora(genes, min_genes = 1)
   
    if (!is.null(rv$ora_result)) {
      shiny::showNotification("ORA complete!", type = "message")
    }
  })
 
  # ORA table
  output$ora_table <- DT::renderDT({
    shiny::req(rv$ora_result)
    df <- rv$ora_result[, c("pathway", "n_overlap", "fold_enrichment", "p_value", "p_adjust")]
    df$fold_enrichment <- round(df$fold_enrichment, 2)
    df$p_value <- format.pval(df$p_value, digits = 2)
    df$p_adjust <- format.pval(df$p_adjust, digits = 2)
    DT::datatable(df, options = list(pageLength = 10))
  })
 
  # ORA plot
  output$plot_ora <- shiny::renderPlot({
    shiny::req(rv$ora_result)
    plot_enrichment_bar(rv$ora_result, show_all = TRUE)
  })
 
  # ML training
  shiny::observeEvent(input$run_ml, {
    shiny::req(rv$scores, rv$clinical, input$ml_outcome)
   
    outcome <- rv$clinical[[input$ml_outcome]]
   
    shiny::withProgress(message = "Training model...", {
      tryCatch({
        rv$ml_model <- death_build_model(
          rv$scores,
          outcome = outcome,
          method = input$ml_method,
          train_ratio = input$ml_train_ratio
        )
        shiny::showNotification("Model training complete!", type = "message")
      }, error = function(e) {
        shiny::showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })
 
  # ML results
  output$ml_results <- shiny::renderPrint({
    shiny::req(rv$ml_model)
    print(rv$ml_model)
  })
 
  # Feature importance plot
  output$plot_importance <- shiny::renderPlot({
    shiny::req(rv$ml_model)
    plots <- plot_model(rv$ml_model, type = "importance")
    print(plots)
  })
 
  # Downloads
  output$dl_scores_csv <- shiny::downloadHandler(
    filename = function() {
      paste0("death_scores_", Sys.Date(), ".csv")
    },
    content = function(file) {
      shiny::req(rv$scores)
      write.csv(rv$scores, file)
    }
  )
 
  output$dl_scores_rds <- shiny::downloadHandler(
    filename = function() {
      paste0("death_scores_", Sys.Date(), ".rds")
    },
    content = function(file) {
      shiny::req(rv$scores)
      saveRDS(rv$scores, file)
    }
  )
 
  output$dl_heatmap <- shiny::downloadHandler(
    filename = function() { "heatmap.pdf" },
    content = function(file) {
      pdf(file, width = 10, height = 8)
      plot_death_heatmap(rv$scores)
      dev.off()
    }
  )
 
  output$dl_boxplot <- shiny::downloadHandler(
    filename = function() { "boxplot.pdf" },
    content = function(file) {
      p <- plot_death_boxplot(rv$scores, group = get_group(),
                              pathways = input$boxplot_pathways)
      ggplot2::ggsave(file, p, width = 10, height = 8)
    }
  )
 
  output$dl_radar <- shiny::downloadHandler(
    filename = function() { "radar.pdf" },
    content = function(file) {
      p <- plot_death_radar(rv$scores, group = get_group())
      ggplot2::ggsave(file, p, width = 8, height = 8)
    }
  )
}


#' Create Minimal Shiny App for Quick Analysis
#'
#' A simplified version of the Shiny app for quick pathway scoring.
#'
#' @export
#'
launch_death_app_mini <- function() {
 
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required")
  }
 
  ui <- shiny::fluidPage(
    shiny::titlePanel("CellDeathAnalysis - Quick Scoring"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::fileInput("file", "Upload Expression (CSV/RDS)"),
        shiny::selectInput("method", "Method",
                           c("zscore", "mean", "median")),
        shiny::actionButton("run", "Calculate", class = "btn-primary"),
        shiny::hr(),
        shiny::downloadButton("download", "Download Results")
      ),
      shiny::mainPanel(
        shiny::verbatimTextOutput("summary"),
        shiny::plotOutput("plot")
      )
    )
  )
 
  server <- function(input, output, session) {
    scores <- shiny::reactiveVal(NULL)
   
    shiny::observeEvent(input$run, {
      shiny::req(input$file)
      ext <- tools::file_ext(input$file$datapath)
      expr <- if (ext == "rds") readRDS(input$file$datapath) else
        as.matrix(read.csv(input$file$datapath, row.names = 1))
     
      scores(calculate_death_score(expr, method = input$method))
    })
   
    output$summary <- shiny::renderPrint({
      shiny::req(scores())
      print(summary(scores()))
    })
   
    output$plot <- shiny::renderPlot({
      shiny::req(scores())
      plot_death_boxplot(scores())
    })
   
    output$download <- shiny::downloadHandler(
      filename = function() "scores.csv",
      content = function(file) write.csv(scores(), file)
    )
  }
 
  shiny::shinyApp(ui, server)
}
