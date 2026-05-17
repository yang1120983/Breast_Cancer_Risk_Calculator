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

      .title-wrapper { padding: 28px 10px 30px 10px; }

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
        box-shadow: 0 4px 12px rgba(0,0,0,0.04),
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

      .gene-value-title {
        font-size: 22px;
        font-weight: 750;
        margin-top: 18px;
        margin-bottom: 8px;
      }

      .abnormal-grid {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 18px;
        margin-top: 20px;
        margin-bottom: 20px;
      }

      .abnormal-card {
        background: #f5f5f7;
        border-radius: 22px;
        padding: 20px;
      }

      .abnormal-title-high {
        color: #ff3b30;
        font-size: 21px;
        font-weight: 800;
        margin-bottom: 12px;
      }

      .abnormal-title-low {
        color: #007aff;
        font-size: 21px;
        font-weight: 800;
        margin-bottom: 12px;
      }

      .abnormal-item {
        display: flex;
        justify-content: space-between;
        padding: 8px 0;
        border-bottom: 1px solid #e5e5ea;
        font-size: 15px;
      }

      .abnormal-item:last-child { border-bottom: none; }

      .abnormal-gene {
        font-weight: 750;
        color: #1d1d1f;
      }

      .abnormal-high {
        color: #ff3b30;
        font-weight: 800;
      }

      .abnormal-low {
        color: #007aff;
        font-weight: 800;
      }

      .dataTables_wrapper { font-size: 14px; }

      table.dataTable thead th {
        background-color: #f5f5f7;
        color: #1d1d1f;
        font-weight: 750;
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
      
      actionButton("load_demo", "Reload Demo Dataset"),
      
      hr(),
      
      uiOutput("sample_selector"),
      
      selectInput(
        "plot_gene_count",
        "Number of genes shown in plot",
        choices = c(
          "10 genes" = 10,
          "20 genes" = 20
        ),
        selected = 10
      ),
      
      actionButton("run_prediction", "Run Prediction", class = "btn-primary"),
      
      hr(),
      
      h4("Gene Check"),
      verbatimTextOutput("gene_check"),
      
      hr(),
      
      uiOutput("download_ui")
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
          div(class = "sample-card-header",
              div(
                h3("Selected Sample Gene Expression"),
                div(class = "sample-subtitle",
                    "Genes shown include the highest- and lowest-expression genes for the selected patient. Genes are ordered from highest to lowest patient expression. Horizontal boxplots show cohort-level distribution, and red triangles indicate the selected patient.")
              )
          ),
          uiOutput("gene_expression_plot_ui"),
          uiOutput("abnormal_gene_summary"),
          div(class = "gene-value-title", "Selected Patient Gene Values"),
          DTOutput("gene_value_table")
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
  
  make_demo_data <- function() {
    set.seed(3888)
    
    data.frame(
      sample_id = paste0("Patient_", 1:100),
      matrix(
        rnorm(100 * 50),
        nrow = 100,
        ncol = 50,
        dimnames = list(NULL, final_top_genes)
      )
    )
  }
  
  dataset <- reactiveVal(make_demo_data())
  prediction_values <- reactiveVal(NULL)
  
  observeEvent(input$load_demo, {
    dataset(make_demo_data())
    prediction_values(NULL)
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
    prediction_values(NULL)
  })
  
  output$sample_selector <- renderUI({
    req(dataset())
    data <- dataset()
    
    selectInput(
      "selected_sample",
      "Select sample / patient",
      choices = data[[1]],
      selected = data[[1]][1]
    )
  })
  
  gene_summary <- reactive({
    req(dataset())
    req(input$selected_sample)
    req(input$plot_gene_count)
    
    data <- dataset()
    
    gene_data <- data[, -1, drop = FALSE]
    gene_data <- as.data.frame(lapply(gene_data, function(x) suppressWarnings(as.numeric(x))))
    
    available_genes <- intersect(final_top_genes, colnames(gene_data))
    
    validate(
      need(length(available_genes) > 1, "Not enough matching genes to create the plot.")
    )
    
    gene_data <- gene_data[, available_genes, drop = FALSE]
    
    selected_index <- which(data[[1]] == input$selected_sample)[1]
    
    all_patient_values <- as.numeric(gene_data[selected_index, available_genes])
    names(all_patient_values) <- available_genes
    
    gene_count <- as.numeric(input$plot_gene_count)
    
    ordered_high <- names(sort(all_patient_values, decreasing = TRUE, na.last = NA))
    ordered_low  <- names(sort(all_patient_values, decreasing = FALSE, na.last = NA))
    
    n_high <- ceiling(gene_count / 2)
    n_low  <- floor(gene_count / 2)
    
    high_genes <- ordered_high[1:min(n_high, length(ordered_high))]
    low_genes  <- ordered_low[1:min(n_low, length(ordered_low))]
    
    top_genes <- unique(c(high_genes, low_genes))
    
    selected_values <- as.numeric(gene_data[selected_index, top_genes])
    
    ordering <- order(selected_values, decreasing = TRUE)
    
    top_genes <- top_genes[ordering]
    selected_values <- selected_values[ordering]
    
    cohort_median <- apply(gene_data[, top_genes, drop = FALSE], 2, median, na.rm = TRUE)
    
    percentile <- sapply(seq_along(top_genes), function(i) {
      gene_values <- gene_data[[top_genes[i]]]
      round(mean(gene_values <= selected_values[i], na.rm = TRUE) * 100, 1)
    })
    
    status <- sapply(percentile, function(p) {
      if (p >= 75) {
        "High"
      } else if (p <= 25) {
        "Low"
      } else {
        "Normal"
      }
    })
    
    list(
      gene_data = gene_data,
      top_genes = top_genes,
      selected_values = selected_values,
      cohort_median = cohort_median,
      percentile = percentile,
      status = status
    )
  })
  
  report_gene_tables <- reactive({
    gs <- gene_summary()
    
    gene_table <- data.frame(
      Gene = gs$top_genes,
      `Patient Value` = round(gs$selected_values, 4),
      `Cohort Median` = round(gs$cohort_median, 4),
      `Patient Percentile` = paste0(gs$percentile, "%"),
      `Expression Status` = gs$status,
      check.names = FALSE
    )
    
    summary_df <- data.frame(
      Gene = gs$top_genes,
      `Patient Value` = round(gs$selected_values, 4),
      Percentile = paste0(gs$percentile, "%"),
      Status = gs$status,
      PercentileNumeric = gs$percentile,
      check.names = FALSE
    )
    
    high_genes <- summary_df[summary_df$Status == "High", ]
    high_genes <- high_genes[order(-high_genes$PercentileNumeric), ]
    high_genes <- head(high_genes[, c("Gene", "Patient Value", "Percentile")], 3)
    
    low_genes <- summary_df[summary_df$Status == "Low", ]
    low_genes <- low_genes[order(low_genes$PercentileNumeric), ]
    low_genes <- head(low_genes[, c("Gene", "Patient Value", "Percentile")], 3)
    
    list(
      gene_table = gene_table,
      high_genes = high_genes,
      low_genes = low_genes
    )
  })
  
  output$gene_expression_plot_ui <- renderUI({
    gene_count <- as.numeric(input$plot_gene_count)
    plot_height <- max(430, 240 + gene_count * 26)
    plotOutput("gene_expression_plot", height = paste0(plot_height, "px"))
  })
  
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
      
      polygon(
        c(outer_r * cos(theta), inner_r * cos(rev(theta))),
        c(outer_r * sin(theta), inner_r * sin(rev(theta))),
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
  
  output$gene_expression_plot <- renderPlot({
    gs <- gene_summary()
    
    top_genes <- gs$top_genes
    selected_values <- gs$selected_values
    
    plot_data <- gs$gene_data[, rev(top_genes), drop = FALSE]
    selected_values <- rev(selected_values)
    
    par(mar = c(5, 8, 4, 2))
    
    boxplot(
      plot_data,
      horizontal = TRUE,
      outline = FALSE,
      col = "#f5f5f7",
      border = "#8e8e93",
      lwd = 1.6,
      las = 1,
      xlab = "Gene Expression",
      main = "Patient Gene Expression Ordered from High to Low",
      cex.axis = 0.9,
      cex.lab = 1.1,
      cex.main = 1.2
    )
    
    points(
      x = selected_values,
      y = seq_along(top_genes),
      pch = 17,
      col = "#ff3b30",
      cex = 1.9
    )
    
    legend(
      "bottomright",
      legend = c("Selected patient"),
      pch = 17,
      col = "#ff3b30",
      bty = "n",
      cex = 0.9
    )
  })
  
  output$abnormal_gene_summary <- renderUI({
    gs <- gene_summary()
    
    summary_df <- data.frame(
      Gene = gs$top_genes,
      PatientValue = gs$selected_values,
      Percentile = gs$percentile,
      Status = gs$status
    )
    
    high_genes <- summary_df[summary_df$Status == "High", ]
    high_genes <- high_genes[order(-high_genes$PatientValue), ]
    high_genes <- head(high_genes, 3)
    
    low_genes <- summary_df[summary_df$Status == "Low", ]
    low_genes <- low_genes[order(low_genes$PatientValue), ]
    low_genes <- head(low_genes, 3)
    
    high_items <- if (nrow(high_genes) == 0) {
      list(div(class = "abnormal-item", span("No high-expression genes detected.")))
    } else {
      lapply(seq_len(nrow(high_genes)), function(i) {
        div(
          class = "abnormal-item",
          span(class = "abnormal-gene", high_genes$Gene[i]),
          span(
            class = "abnormal-high",
            paste0(
              "▲ ",
              round(high_genes$PatientValue[i], 4),
              " (",
              high_genes$Percentile[i],
              "%)"
            )
          )
        )
      })
    }
    
    low_items <- if (nrow(low_genes) == 0) {
      list(div(class = "abnormal-item", span("No low-expression genes detected.")))
    } else {
      lapply(seq_len(nrow(low_genes)), function(i) {
        div(
          class = "abnormal-item",
          span(class = "abnormal-gene", low_genes$Gene[i]),
          span(
            class = "abnormal-low",
            paste0(
              "▼ ",
              round(low_genes$PatientValue[i], 4),
              " (",
              low_genes$Percentile[i],
              "%)"
            )
          )
        )
      })
    }
    
    div(
      class = "abnormal-grid",
      div(
        class = "abnormal-card",
        div(class = "abnormal-title-high", "Top High Expression Genes"),
        high_items
      ),
      div(
        class = "abnormal-card",
        div(class = "abnormal-title-low", "Top Low Expression Genes"),
        low_items
      )
    )
  })
  
  output$gene_value_table <- renderDT({
    gs <- gene_summary()
    
    expression_status <- sapply(gs$status, function(s) {
      if (s == "High") {
        "<span style='color:#ff3b30; font-weight:800;'>▲ High</span>"
      } else if (s == "Low") {
        "<span style='color:#007aff; font-weight:800;'>▼ Low</span>"
      } else {
        "<span style='color:#8e8e93; font-weight:800;'>● Normal</span>"
      }
    })
    
    value_table <- data.frame(
      Gene = gs$top_genes,
      `Patient Value` = round(gs$selected_values, 4),
      `Cohort Median` = round(gs$cohort_median, 4),
      `Patient Percentile` = paste0(gs$percentile, "%"),
      `Expression Status` = expression_status,
      check.names = FALSE
    )
    
    dt <- datatable(
      value_table,
      escape = FALSE,
      rownames = FALSE,
      options = list(
        pageLength = min(as.numeric(input$plot_gene_count), 10),
        lengthMenu = c(10, 20),
        dom = "tip",
        ordering = FALSE,
        autoWidth = TRUE
      )
    )
    
    dt <- formatStyle(dt, "Gene", fontWeight = "700")
    dt <- formatStyle(dt, "Patient Value", fontWeight = "700")
    
    dt
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
    
    patient_seed <- sum(utf8ToInt(input$selected_sample))
    set.seed(patient_seed)
    
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
      div(
        class = low_class,
        div(class = "risk-level-title", "Low"),
        div(class = "risk-level-range", "0.0 – 0.4")
      ),
      div(
        class = medium_class,
        div(class = "risk-level-title", "Medium"),
        div(class = "risk-level-range", "0.4 – 0.7")
      ),
      div(
        class = high_class,
        div(class = "risk-level-title", "High"),
        div(class = "risk-level-range", "0.7 – 1.0")
      )
    )
  })
  
  output$uncertainty_box <- renderUI({
    req(prediction_values())
    
    div(
      class = "uncertainty-pill",
      paste0(prediction_values()$ci_low, " - ", prediction_values()$ci_high)
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
          style = "color:#ff3b30; font-weight:850; font-size:52px; margin-top:0px; margin-bottom:18px;"
        ),
        p(
          "The gene expression profile indicates a high predicted risk.",
          style = "font-size:22px; font-weight:500; margin-bottom:12px; line-height:1.45;"
        ),
        p(
          "Further clinical assessment and molecular subtype assessment is recommended.",
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
  
  output$download_ui <- renderUI({
    req(prediction_values())
    
    downloadButton(
      "download_report",
      "Download Patient Report",
      class = "btn-primary"
    )
  })
  
  output$download_report <- downloadHandler(
    filename = function() {
      paste0(input$selected_sample, "_risk_report_", Sys.Date(), ".pdf")
    },
    
    content = function(file) {
      req(prediction_values())
      req(input$selected_sample)
      req(gene_summary())
      
      vals <- prediction_values()
      report_tables <- report_gene_tables()
      
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
          gene_table = report_tables$gene_table,
          high_genes = report_tables$high_genes,
          low_genes = report_tables$low_genes
        ),
        envir = new.env(parent = globalenv())
      )
      
      file.copy(output_pdf, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui, server)