library(shiny)
library(DT)
library(rmarkdown)

options(shiny.maxRequestSize = 100 * 1024^2)

ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("
      body {
        background-color: #f5f5f7;
        font-family: -apple-system, BlinkMacSystemFont, 'SF Pro Display', 'Helvetica Neue', Arial, sans-serif;
        color: #1d1d1f;
      }

      .title-wrapper {
        padding: 28px 10px 30px 10px;
      }

      .app-title {
        font-size: 56px;
        font-weight: 800;
        letter-spacing: -0.5px;
        color: #111111;
        margin-bottom: 6px;
      }

      .app-subtitle {
        font-size: 25px;
        color: #6e6e73;
        font-weight: 450;
      }

      .card {
        background-color: #ffffff;
        border-radius: 32px;
        padding: 28px;
        box-shadow:
          0 4px 12px rgba(0,0,0,0.04),
          0 12px 32px rgba(0,0,0,0.06);
        margin-bottom: 22px;
        border: none;
      }

      .metric-card {
        height: 320px;
        display: flex;
        flex-direction: column;
      }

      .metric-card h4 {
        margin-top: 0;
        font-size: 24px;
        font-weight: 750;
        color: #1d1d1f;
      }

      .metric-content {
        flex: 1;
        display: flex;
        flex-direction: column;
        justify-content: center;
      }

      h3, h4 {
        font-weight: 700;
        color: #1d1d1f;
      }

      .uncertainty-pill {
        background: rgba(0,122,255,0.12);
        color: #007aff;
        border-radius: 999px;
        padding: 22px 26px;
        font-size: 32px;
        font-weight: 800;
        text-align: center;
        letter-spacing: -1.5px;
      }

      .uncertainty-note {
        text-align: center;
        margin-top: 14px;
        color: #6e6e73;
        font-size: 15px;
        font-weight: 500;
      }

      .warning-box {
        background-color: #fff9e6;
        border-radius: 20px;
        padding: 18px;
        color: #6e4b00;
        border: none;
      }

      .btn-primary {
        background-color: #007aff !important;
        border-color: #007aff !important;
        border-radius: 14px !important;
        font-weight: 650;
      }

      .form-control {
        border-radius: 14px;
        border: 1px solid #d2d2d7;
      }

      .well {
        background-color: #ffffff;
        border-radius: 24px;
        border: none;
        box-shadow: 0 4px 18px rgba(0,0,0,0.05);
      }

      .risk-level-list {
        display: flex;
        flex-direction: column;
        gap: 10px;
      }

      .risk-level-card {
        border-radius: 18px;
        padding: 12px 16px;
        display: flex;
        justify-content: space-between;
        align-items: center;
        opacity: 0.55;
      }

      .risk-level-active {
        opacity: 1;
        transform: scale(1.02);
        box-shadow: 0 8px 18px rgba(0,0,0,0.08);
      }

      .risk-level-title {
        font-size: 20px;
        font-weight: 800;
      }

      .risk-level-range {
        font-size: 14px;
        font-weight: 700;
      }

      .risk-level-low {
        background: rgba(52,199,89,0.14);
        color: #34c759;
      }

      .risk-level-medium {
        background: rgba(255,204,0,0.20);
        color: #9a7a00;
      }

      .risk-level-high {
        background: rgba(255,59,48,0.12);
        color: #ff3b30;
      }

      .model-grid {
        display: grid;
        grid-template-columns: repeat(4, 1fr);
        gap: 18px;
        margin-top: 18px;
      }

      .model-card {
        background: #f5f5f7;
        border-radius: 22px;
        padding: 20px;
        text-align: center;
      }

      .model-name {
        font-size: 16px;
        font-weight: 700;
        color: #1d1d1f;
        margin-bottom: 12px;
      }

      .model-value {
        font-size: 34px;
        font-weight: 850;
        color: #007aff;
        letter-spacing: -1px;
      }

      .model-role {
        font-size: 13px;
        color: #6e6e73;
        margin-top: 10px;
      }

      .sample-card-header {
        display: flex;
        justify-content: space-between;
        align-items: center;
        margin-bottom: 18px;
      }

      .sample-subtitle {
        color: #6e6e73;
        font-size: 15px;
        margin-top: 4px;
      }

      .dataTables_wrapper {
        font-size: 15px;
      }

      table.dataTable {
        border-collapse: collapse !important;
        border-radius: 18px;
        overflow: hidden;
      }

      table.dataTable thead th {
        background-color: #f5f5f7;
        color: #1d1d1f;
        font-weight: 750;
        border-bottom: none !important;
      }

      table.dataTable tbody td {
        border-top: 1px solid #eeeeee;
      }
    "))
  ),
  
  div(class = "title-wrapper",
      div(class = "app-title", "Breast Cancer Risk Calculator"),
      div(class = "app-subtitle", "Clinical Risk Prediction Dashboard")
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      h4("Input"),
      
      fileInput(
        "demo_file",
        "Upload demo dataset CSV",
        accept = c(".csv", ".txt")
      ),
      
      actionButton("load_demo", "Load Demo Dataset"),
      
      hr(),
      
      uiOutput("sample_selector"),
      
      selectInput(
        "gene_panel",
        "Gene panel size",
        choices = c(
          "50 genes" = 50,
          "100 genes" = 100,
          "200 genes" = 200,
          "500 genes" = 500,
          "700 genes" = 700
        ),
        selected = 50
      ),
      
      actionButton("run_prediction", "Run Prediction", class = "btn-primary"),
      
      hr(),
      
      h4("Gene Check"),
      verbatimTextOutput("gene_check"),
      
      hr(),
      
      downloadButton(
        "download_report",
        "Download Patient Report",
        class = "btn-primary"
      )
    ),
    
    mainPanel(
      
      fluidRow(
        column(
          4,
          div(class = "card metric-card",
              h4("Main Risk Score"),
              div(class = "metric-content",
                  plotOutput("risk_gauge", height = "250px")
              )
          )
        ),
        
        column(
          4,
          div(class = "card metric-card",
              h4("Risk Category"),
              div(class = "metric-content",
                  uiOutput("risk_category_levels")
              )
          )
        ),
        
        column(
          4,
          div(class = "card metric-card",
              h4("Uncertainty"),
              div(class = "metric-content",
                  uiOutput("uncertainty_box"),
                  div(class = "uncertainty-note", "95% confidence interval")
              )
          )
        )
      ),
      
      div(class = "card",
          h3("Risk Interpretation"),
          uiOutput("risk_interpretation_box")
      ),
      
      div(class = "card",
          h3("Base Model Predictions"),
          uiOutput("base_model_cards")
      ),
      
      div(class = "card",
          div(class = "sample-card-header",
              div(
                h3("Selected Sample Gene Expression"),
                div(class = "sample-subtitle",
                    "Expression values are rounded to 4 decimal places for readability.")
              )
          ),
          DTOutput("sample_table")
      ),
      
      div(class = "card",
          h3("Clinical Warning / Limitation"),
          div(class = "warning-box", verbatimTextOutput("warning_box"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  final_top_genes <- paste0("Gene", 1:50)
  
  dataset <- reactiveVal(NULL)
  prediction_values <- reactiveVal(NULL)
  
  output$risk_gauge <- renderPlot({
    req(prediction_values())
    
    risk <- prediction_values()$meta
    
    par(mar = c(0, 0, 0, 0))
    
    plot(
      0, 0,
      type = "n",
      xlim = c(-1.15, 1.15),
      ylim = c(-0.45, 1.15),
      axes = FALSE,
      xlab = "",
      ylab = "",
      asp = 1
    )
    
    draw_arc <- function(start, end, col) {
      theta <- seq(start, end, length.out = 300)
      
      outer_r <- 1
      inner_r <- 0.6
      
      x_outer <- outer_r * cos(theta)
      y_outer <- outer_r * sin(theta)
      x_inner <- inner_r * cos(rev(theta))
      y_inner <- inner_r * sin(rev(theta))
      
      polygon(
        c(x_outer, x_inner),
        c(y_outer, y_inner),
        col = col,
        border = NA
      )
    }
    
    draw_arc(pi, pi * 0.60, "#34c759")
    draw_arc(pi * 0.60, pi * 0.30, "#ffcc00")
    draw_arc(pi * 0.30, 0, "#ff3b30")
    
    angle <- pi * (1 - risk)
    needle_length <- 0.78
    
    segments(
      0, 0,
      needle_length * cos(angle),
      needle_length * sin(angle),
      lwd = 5,
      col = "#1d1d1f"
    )
    
    points(0, 0, pch = 19, cex = 2.2, col = "#1d1d1f")
    
    text(-1.03, -0.08, "0", cex = 1.15, font = 2, col = "#1d1d1f")
    text(1.03, -0.08, "1", cex = 1.15, font = 2, col = "#1d1d1f")
    
    text(0, -0.27, risk, cex = 2.5, font = 2, col = "#ff3b30")
    text(0, -0.47, "Risk Score", cex = 0.95, font = 2, col = "#6e6e73")
    text(
      0,
      -0.6,
      "Range: 0 (Low Risk) – 1 (High Risk)",
      cex = 0.95,
      font = 2,
      col = "#6e6e73"
    )
  })
  
  observeEvent(input$demo_file, {
    req(input$demo_file)
    
    data <- read.csv(input$demo_file$datapath, check.names = FALSE)
    
    data$PatientID <- paste0(
      "Patient_",
      sprintf("%03d", 1:nrow(data))
    )
    
    data <- data[, c(ncol(data), 1:(ncol(data) - 1))]
    dataset(data)
  })
  
  observeEvent(input$load_demo, {
    
    set.seed(3888)
    
    demo_data <- data.frame(
      sample_id = paste0("Patient_", 1:100),
      matrix(
        rnorm(100 * 50),
        nrow = 100,
        ncol = 50,
        dimnames = list(NULL, final_top_genes)
      )
    )
    
    dataset(demo_data)
  })
  
  output$sample_selector <- renderUI({
    req(dataset())
    data <- dataset()
    
    selectInput(
      "selected_sample",
      "Select sample / patient",
      choices = data[[1]]
    )
  })
  
  output$gene_check <- renderText({
    req(dataset())
    
    data <- dataset()
    available_genes <- colnames(data)
    
    found_genes <- intersect(final_top_genes, available_genes)
    missing_genes <- setdiff(final_top_genes, available_genes)
    
    paste(
      "Required genes:", length(final_top_genes),
      "\nFound genes:", length(found_genes),
      "\nMissing genes:", length(missing_genes),
      "\n\nMissing gene list:",
      ifelse(length(missing_genes) == 0, "None", paste(missing_genes, collapse = ", "))
    )
  })
  
  observeEvent(input$run_prediction, {
    req(dataset())
    req(input$selected_sample)
    
    rf_pred <- round(runif(1, 0.65, 0.90), 3)
    ridge_pred <- round(runif(1, 0.60, 0.88), 3)
    nb_pred <- round(runif(1, 0.58, 0.86), 3)
    
    meta_risk <- round(0.45 * rf_pred + 0.35 * ridge_pred + 0.20 * nb_pred, 3)
    
    ci_low <- max(0, round(meta_risk - runif(1, 0.05, 0.10), 3))
    ci_high <- min(1, round(meta_risk + runif(1, 0.05, 0.10), 3))
    
    prediction_values(list(
      rf = rf_pred,
      ridge = ridge_pred,
      nb = nb_pred,
      meta = meta_risk,
      ci_low = ci_low,
      ci_high = ci_high
    ))
  })
  
  output$risk_category_levels <- renderUI({
    req(prediction_values())
    
    risk <- prediction_values()$meta
    
    low_class <- if (risk < 0.4) {
      "risk-level-card risk-level-low risk-level-active"
    } else {
      "risk-level-card risk-level-low"
    }
    
    medium_class <- if (risk >= 0.4 && risk < 0.7) {
      "risk-level-card risk-level-medium risk-level-active"
    } else {
      "risk-level-card risk-level-medium"
    }
    
    high_class <- if (risk >= 0.7) {
      "risk-level-card risk-level-high risk-level-active"
    } else {
      "risk-level-card risk-level-high"
    }
    
    div(
      class = "risk-level-list",
      
      div(class = low_class,
          div(class = "risk-level-title", "Low"),
          div(class = "risk-level-range", "0.0 – 0.4")
      ),
      
      div(class = medium_class,
          div(class = "risk-level-title", "Medium"),
          div(class = "risk-level-range", "0.4 – 0.7")
      ),
      
      div(class = high_class,
          div(class = "risk-level-title", "High"),
          div(class = "risk-level-range", "0.7 – 1.0")
      )
    )
  })
  
  output$uncertainty_box <- renderUI({
    req(prediction_values())
    
    div(
      class = "uncertainty-pill",
      paste0(
        prediction_values()$ci_low,
        " - ",
        prediction_values()$ci_high
      )
    )
  })
  
  output$risk_interpretation_box <- renderUI({
    req(prediction_values())
    
    risk <- prediction_values()$meta
    
    if (risk >= 0.7) {
      div(
        style = "
          background: rgba(255,59,48,0.08);
          border: 1px solid rgba(255,59,48,0.20);
          border-radius: 22px;
          padding: 18px 26px 20px 26px;
        ",
        h2(
          "⚠️ HIGH RISK",
          style = "
            color:#ff3b30;
            font-weight:850;
            font-size:52px;
            margin-top:0px;
            margin-bottom:18px;
          "
        ),
        p(
          "The gene expression profile indicates a high predicted risk.",
          style = "font-size:22px; font-weight:500; margin-bottom:12px; line-height:1.45;"
        ),
        p(
          "Further clinical assessment and molecular subtyping are recommended.",
          style = "font-size:22px; font-weight:500; line-height:1.45; margin-bottom:0px;"
        )
      )
    } else if (risk >= 0.4) {
      div(
        style = "
          background: rgba(255,149,0,0.10);
          border: 1px solid rgba(255,149,0,0.25);
          border-radius: 22px;
          padding: 18px 26px 20px 26px;
        ",
        h2(
          "⚠️ MEDIUM RISK",
          style = "color:#ff9500; font-weight:850; font-size:52px; margin-top:0px; margin-bottom:18px;"
        ),
        p(
          "The gene expression profile indicates an intermediate predicted risk.",
          style = "font-size:22px; font-weight:500; margin-bottom:12px; line-height:1.45;"
        ),
        p(
          "Additional clinical review may be useful before making decisions.",
          style = "font-size:22px; font-weight:500; line-height:1.45; margin-bottom:0px;"
        )
      )
    } else {
      div(
        style = "
          background: rgba(52,199,89,0.10);
          border: 1px solid rgba(52,199,89,0.25);
          border-radius: 22px;
          padding: 18px 26px 20px 26px;
        ",
        h2(
          "✅ LOW RISK",
          style = "color:#34c759; font-weight:850; font-size:52px; margin-top:0px; margin-bottom:18px;"
        ),
        p(
          "The gene expression profile indicates a low predicted risk.",
          style = "font-size:22px; font-weight:500; margin-bottom:12px; line-height:1.45;"
        ),
        p(
          "Routine monitoring may be appropriate depending on clinical context.",
          style = "font-size:22px; font-weight:500; line-height:1.45; margin-bottom:0px;"
        )
      )
    }
  })
  
  output$base_model_cards <- renderUI({
    req(prediction_values())
    
    vals <- prediction_values()
    
    div(
      class = "model-grid",
      
      div(class = "model-card",
          div(class = "model-name", "Random Forest"),
          div(class = "model-value", vals$rf),
          div(class = "model-role", "Base model")
      ),
      
      div(class = "model-card",
          div(class = "model-name", "Ridge Regression"),
          div(class = "model-value", vals$ridge),
          div(class = "model-role", "Base model")
      ),
      
      div(class = "model-card",
          div(class = "model-name", "Naive Bayes"),
          div(class = "model-value", vals$nb),
          div(class = "model-role", "Base model")
      ),
      
      div(class = "model-card",
          div(class = "model-name", "Attention Meta-model"),
          div(class = "model-value", vals$meta),
          div(class = "model-role", "Final risk integration")
      )
    )
  })
  
  output$sample_table <- renderDT({
    req(dataset())
    req(input$selected_sample)
    
    data <- dataset()
    selected_row <- data[data[[1]] == input$selected_sample, ]
    
    selected_row[-1] <- round(selected_row[-1], 4)
    
    datatable(
      selected_row,
      options = list(
        scrollX = TRUE,
        pageLength = 5
      )
    )
  })
  
  output$warning_box <- renderText({
    paste(
      "This is currently a front-end demonstration module.",
      "The final version will replace placeholder values with trained RF, Ridge, NB, and attention meta-model outputs.",
      "The uploaded dataset must contain the required final top 50 genes.",
      "Prediction uncertainty will be estimated using bootstrap models.",
      "Clinical interpretation should be cautious if genes are missing or if the demo dataset differs from the training datasets.",
      sep = "\n"
    )
  })
  
  output$download_report <- downloadHandler(
    
    filename = function() {
      paste0(
        input$selected_sample,
        "_risk_report_",
        Sys.Date(),
        ".pdf"
      )
    },
    
    content = function(file) {
      
      req(prediction_values())
      req(input$selected_sample)
      
      vals <- prediction_values()
      
      risk_category <- ifelse(
        vals$meta >= 0.7,
        "High Risk",
        ifelse(vals$meta >= 0.4, "Medium Risk", "Low Risk")
      )
      
      temp_report <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", temp_report, overwrite = TRUE)
      
      output_pdf <- rmarkdown::render(
        input = temp_report,
        output_format = "pdf_document",
        output_file = "patient_risk_report.pdf",
        output_dir = tempdir(),
        params = list(
          patient_id = input$selected_sample,
          risk_score = vals$meta,
          risk_category = risk_category,
          ci_low = vals$ci_low,
          ci_high = vals$ci_high,
          rf = vals$rf,
          ridge = vals$ridge,
          nb = vals$nb
        ),
        envir = new.env(parent = globalenv())
      )
      
      file.copy(output_pdf, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui, server)