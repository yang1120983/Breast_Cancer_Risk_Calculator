library(shiny)
library(DT)
library(rmarkdown)
library(GEOquery)
library(Biobase)
library(limma)
library(randomForest)
library(glmnet)
library(nnet)
library(class)

options(shiny.maxRequestSize = 100 * 1024^2)

find_multiclass_model_dir <- function() {
  candidates <- c(
    file.path("Pre-processing", "models", "multiclass"),
    file.path("..", "Pre-processing", "models", "multiclass"),
    file.path("..", "..", "Pre-processing", "models", "multiclass")
  )
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) == 0) return(NULL)
  existing[1]
}

load_multiclass_model <- function() {
  model_dir <- find_multiclass_model_dir()
  if (is.null(model_dir)) {
    return(list(available = FALSE, message = "Could not find Pre-processing/models/multiclass."))
  }
  
  required_files <- c(
    "fit_rf.rds", "fit_lasso.rds", "fit_knn.rds",
    "panel_genes.rds", "imputation_params.rds",
    "scale_params.rds", "batch_correction_params.rds"
  )
  missing_files <- required_files[!file.exists(file.path(model_dir, required_files))]
  if (length(missing_files) > 0) {
    return(list(
      available = FALSE,
      message = paste("Missing multiclass model files:", paste(missing_files, collapse = ", "))
    ))
  }
  
  model_disagreement_file <- file.path(model_dir, "model_disagreement.rds")
  bootstrap_risk_file <- file.path(model_dir, "bootstrap_risk_scores.rds")
  
  list(
    available = TRUE,
    model_dir = normalizePath(model_dir),
    fit_rf = readRDS(file.path(model_dir, "fit_rf.rds")),
    fit_lasso = readRDS(file.path(model_dir, "fit_lasso.rds")),
    fit_knn = readRDS(file.path(model_dir, "fit_knn.rds")),
    fit_multinom = if (file.exists(file.path(model_dir, "fit_multinom.rds"))) {
      readRDS(file.path(model_dir, "fit_multinom.rds"))
    } else {
      NULL
    },
    panel_genes = readRDS(file.path(model_dir, "panel_genes.rds")),
    imputation = readRDS(file.path(model_dir, "imputation_params.rds")),
    scale = readRDS(file.path(model_dir, "scale_params.rds")),
    batch = readRDS(file.path(model_dir, "batch_correction_params.rds")),
    model_disagreement = if (file.exists(model_disagreement_file)) readRDS(model_disagreement_file) else NULL,
    bootstrap_risk_scores = if (file.exists(bootstrap_risk_file)) readRDS(bootstrap_risk_file) else NULL,
    uncertainty_rds_name = paste(
      basename(c(model_disagreement_file, bootstrap_risk_file)[file.exists(c(model_disagreement_file, bootstrap_risk_file))]),
      collapse = ", "
    )
  )
}

find_geo_cache_dir <- function() {
  candidates <- c(
    file.path("Pre-processing", "data", "geo-cache"),
    file.path("..", "Pre-processing", "data", "geo-cache"),
    file.path("..", "..", "Pre-processing", "data", "geo-cache")
  )
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) == 0) return(NULL)
  existing[1]
}

prepare_geoquery_platform_cache <- function(target_dir) {
  cache_dir <- find_geo_cache_dir()
  if (is.null(cache_dir)) return(invisible(FALSE))
  
  gpl_files <- list.files(cache_dir, pattern = "^GPL.*\\.soft\\.gz$", full.names = TRUE)
  for (src in gpl_files) {
    dest <- file.path(target_dir, basename(src))
    if (!file.exists(dest)) {
      linked <- suppressWarnings(file.symlink(normalizePath(src), dest))
      if (!isTRUE(linked)) file.copy(src, dest, overwrite = TRUE)
    }
  }
  
  invisible(TRUE)
}

standardize_metadata_ids <- function(meta) {
  if (is.null(meta) || nrow(meta) == 0) return(meta)
  
  names_lower <- tolower(names(meta))
  sample_candidates <- c("geo_accession", "sample_id", "specimen_id", "gsm", "id")
  sample_col <- names(meta)[match(sample_candidates, names_lower, nomatch = 0)]
  sample_col <- sample_col[nzchar(sample_col)][1]
  
  if (is.na(sample_col) || is.null(sample_col)) {
    meta$sample_id <- rownames(meta)
  } else if (!"sample_id" %in% names(meta)) {
    meta$sample_id <- as.character(meta[[sample_col]])
  }
  
  if (!"sample_id" %in% names(meta) || any(is.na(meta$sample_id) | meta$sample_id == "")) {
    meta$sample_id <- sprintf("SAMPLE%03d", seq_len(nrow(meta)))
  }
  
  meta$patient_display_id <- sprintf("Patient %d", seq_len(nrow(meta)))
  meta
}

first_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  vapply(x, function(value) {
    if (is.na(value) || !nzchar(value)) return(NA_character_)
    parts <- unlist(strsplit(value, "///|//|;|,|\\|", perl = TRUE))
    trimws(parts[1])
  }, character(1))
}

collapse_expression_to_genes <- function(eset) {
  expr <- Biobase::exprs(eset)
  fdat <- Biobase::fData(eset)
  
  expr <- apply(expr, 2, function(x) suppressWarnings(as.numeric(x)))
  rownames(expr) <- rownames(Biobase::exprs(eset))
  
  if (stats::quantile(expr, 0.99, na.rm = TRUE) > 100) {
    expr <- log2(expr + 1)
  }
  
  expr <- limma::normalizeBetweenArrays(expr, method = "quantile")
  
  symbol_col <- grep("gene.?symbol|^symbol$", names(fdat), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(symbol_col) || is.null(symbol_col)) {
    gene_symbols <- rownames(expr)
  } else {
    gene_symbols <- first_gene_symbol(fdat[[symbol_col]])
  }
  
  missing_symbol <- is.na(gene_symbols) | !nzchar(gene_symbols)
  gene_symbols[missing_symbol] <- rownames(expr)[missing_symbol]
  
  limma::avereps(expr, ID = gene_symbols)
}

read_geo_series_matrix_upload <- function(file_info) {
  safe_name <- gsub("[^A-Za-z0-9_.-]", "_", basename(file_info$name))
  upload_copy <- file.path(tempdir(), safe_name)
  file.copy(file_info$datapath, upload_copy, overwrite = TRUE)
  prepare_geoquery_platform_cache(dirname(upload_copy))
  
  geo <- GEOquery::getGEO(filename = upload_copy, GSEMatrix = TRUE)
  if (is.list(geo)) geo <- geo[[1]]
  if (!inherits(geo, "ExpressionSet")) {
    stop("The uploaded file could not be read as a GEO series-matrix ExpressionSet.")
  }
  
  meta <- standardize_metadata_ids(as.data.frame(Biobase::pData(geo), stringsAsFactors = FALSE))
  expr_gene <- collapse_expression_to_genes(geo)
  
  expr_table <- as.data.frame(t(expr_gene), check.names = FALSE)
  expr_table <- data.frame(
    sample_id = rownames(expr_table),
    expr_table,
    check.names = FALSE
  )
  
  list(
    data = expr_table,
    source = "GEO series-matrix TXT.GZ",
    n_probes = nrow(Biobase::exprs(geo)),
    n_genes = nrow(expr_gene),
    n_samples = ncol(expr_gene),
    platform = paste(unique(Biobase::annotation(geo)), collapse = ", ")
  )
}

read_tabular_upload <- function(file_info) {
  file_name <- tolower(file_info$name)
  connection <- if (grepl("\\.gz$", file_name)) gzfile(file_info$datapath, open = "rt") else file_info$datapath
  on.exit(if (inherits(connection, "connection")) close(connection), add = TRUE)
  
  if (grepl("\\.csv(\\.gz)?$", file_name)) {
    data <- read.csv(connection, check.names = FALSE)
  } else {
    data <- read.delim(connection, check.names = FALSE, comment.char = "", quote = "\"")
  }
  
  if (!"PatientID" %in% names(data)) {
    data$PatientID <- paste0("Patient_", sprintf("%03d", seq_len(nrow(data))))
    data <- data[, c(ncol(data), seq_len(ncol(data) - 1)), drop = FALSE]
  }
  
  list(
    data = data,
    source = "tabular upload",
    n_samples = nrow(data),
    n_genes = max(0, ncol(data) - 1),
    platform = "Not available"
  )
}

read_expression_upload <- function(file_info) {
  file_name <- tolower(file_info$name)
  looks_like_geo <- grepl("\\.txt\\.gz$|series.*matrix|gse", file_name)
  
  if (looks_like_geo) {
    geo_result <- tryCatch(
      read_geo_series_matrix_upload(file_info),
      error = function(e) {
        warning("GEO parser failed, trying as a regular table: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(geo_result)) return(geo_result)
  }
  
  read_tabular_upload(file_info)
}

scale_with_training_params <- function(x, scale_params) {
  mu <- scale_params$mu[names(x)]
  sd <- scale_params$sd[names(x)]
  sd[is.na(sd) | sd == 0] <- 1
  x_scaled <- sweep(as.matrix(x), 2, mu, "-")
  x_scaled <- sweep(x_scaled, 2, sd, "/")
  as.data.frame(x_scaled)
}

align_probability_vector <- function(prob_vec, subtype_levels) {
  out <- setNames(rep(NA_real_, length(subtype_levels)), subtype_levels)
  common <- intersect(names(prob_vec), subtype_levels)
  out[common] <- prob_vec[common]
  out
}

predict_multiclass_model_probabilities <- function(x_i, x_s, multiclass_model, subtype_levels) {
  probs <- list()
  
  rf_prob <- predict(multiclass_model$fit_rf, x_i, type = "prob")
  probs$RF <- align_probability_vector(
    stats::setNames(as.numeric(rf_prob[1, ]), colnames(rf_prob)),
    subtype_levels
  )
  
  lasso_prob <- tryCatch({
    pr <- predict(multiclass_model$fit_lasso, as.matrix(x_s), type = "response", s = "lambda.min")
    if (length(dim(pr)) == 3) pr <- pr[, , 1]
    if (is.null(dim(pr))) pr <- t(as.matrix(pr))
    align_probability_vector(stats::setNames(as.numeric(pr[1, ]), colnames(pr)), subtype_levels)
  }, error = function(e) NULL)
  if (!is.null(lasso_prob)) probs$LassoMC <- lasso_prob
  
  if (!is.null(multiclass_model$fit_multinom)) {
    multinom_prob <- tryCatch({
      pm <- predict(multiclass_model$fit_multinom, newdata = data.frame(x_s), type = "probs")
      if (is.null(dim(pm))) pm <- t(as.matrix(pm))
      align_probability_vector(stats::setNames(as.numeric(pm[1, ]), colnames(pm)), subtype_levels)
    }, error = function(e) NULL)
    if (!is.null(multinom_prob)) probs$Multinom <- multinom_prob
  }
  
  knn_prob <- tryCatch({
    fit <- multiclass_model$fit_knn
    common <- intersect(colnames(fit$X_train), colnames(x_s))
    pred <- class::knn(
      train = fit$X_train[, common, drop = FALSE],
      test = as.matrix(x_s[, common, drop = FALSE]),
      cl = fit$y_train,
      k = fit$k,
      prob = TRUE
    )
    winning_prob <- attr(pred, "prob")
    winning_class <- as.character(pred)
    levels_knn <- levels(fit$y_train)
    pm <- setNames(rep((1 - winning_prob) / max(1, length(levels_knn) - 1), length(levels_knn)), levels_knn)
    pm[winning_class] <- winning_prob
    align_probability_vector(pm, subtype_levels)
  }, error = function(e) NULL)
  if (!is.null(knn_prob)) probs$kNN <- knn_prob
  
  probs
}

calculate_multiclass_uncertainty <- function(model_probabilities, model_disagreement_reference = NULL) {
  prob_stack <- do.call(rbind, model_probabilities)
  avg_probs <- colMeans(prob_stack, na.rm = TRUE)
  avg_probs <- avg_probs / sum(avg_probs, na.rm = TRUE)
  pred_subtype <- names(avg_probs)[which.max(avg_probs)]
  disagreement <- stats::sd(prob_stack[, pred_subtype], na.rm = TRUE)
  entropy <- {
    p <- avg_probs[avg_probs > 0]
    -sum(p * log(p))
  }
  confidence <- max(avg_probs, na.rm = TRUE)
  
  if (!is.null(model_disagreement_reference) && "disagreement" %in% names(model_disagreement_reference)) {
    low_disagreement_cutoff <- stats::quantile(model_disagreement_reference$disagreement, 0.25, na.rm = TRUE)
    high_disagreement_cutoff <- stats::quantile(model_disagreement_reference$disagreement, 0.75, na.rm = TRUE)
  } else {
    low_disagreement_cutoff <- 0.08
    high_disagreement_cutoff <- 0.15
  }
  
  level <- if (confidence >= 0.70 && disagreement <= low_disagreement_cutoff) {
    "Low"
  } else if (confidence >= 0.50 && disagreement <= high_disagreement_cutoff) {
    "Moderate"
  } else {
    "High"
  }
  
  list(
    level = level,
    confidence = confidence,
    disagreement = disagreement,
    entropy = entropy,
    low_disagreement_cutoff = low_disagreement_cutoff,
    high_disagreement_cutoff = high_disagreement_cutoff,
    model_count = nrow(prob_stack),
    model_probabilities = prob_stack,
    average_probabilities = avg_probs
  )
}

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

      .subtype-result-grid {
        display: grid;
        grid-template-columns: 1.1fr 1fr;
        gap: 20px;
        margin-top: 12px;
      }

      .subtype-primary {
        background: #f5f5f7;
        border-radius: 24px;
        padding: 22px;
      }

      .subtype-label {
        color: #6e6e73;
        font-size: 14px;
        font-weight: 750;
        text-transform: uppercase;
        letter-spacing: 0.04em;
      }

      .subtype-name {
        font-size: 42px;
        font-weight: 850;
        letter-spacing: -0.8px;
        color: #1d1d1f;
        margin-top: 4px;
      }

      .subtype-prob {
        font-size: 34px;
        font-weight: 850;
        color: #007aff;
        margin-top: 4px;
      }

      .subtype-pill {
        display: inline-block;
        margin-top: 14px;
        border-radius: 999px;
        padding: 9px 14px;
        background: rgba(0,122,255,0.12);
        color: #007aff;
        font-weight: 800;
      }

      .subtype-prob-row {
        display: grid;
        grid-template-columns: 95px 1fr 58px;
        gap: 12px;
        align-items: center;
        margin: 12px 0;
        font-size: 15px;
        font-weight: 700;
      }

      .subtype-bar-track {
        height: 12px;
        background: #e5e5ea;
        border-radius: 999px;
        overflow: hidden;
      }

      .subtype-bar-fill {
        height: 100%;
        border-radius: 999px;
        background: #007aff;
      }

      .subtype-muted {
        color: #6e6e73;
        font-size: 15px;
        line-height: 1.45;
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
        "Upload expression dataset CSV, TXT, or TXT.GZ",
        accept = c(".csv", ".txt", ".gz", ".txt.gz")
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
      
      uiOutput("multiclass_section"),
      
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
  multiclass_model <- load_multiclass_model()
  
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
  upload_summary <- reactiveVal(NULL)
  prediction_values <- reactiveVal(NULL)
  
  patient_lookup <- reactive({
    req(dataset())
    data <- dataset()
    data.frame(
      patient_display_id = paste0("Patient_", seq_len(nrow(data))),
      sample_id = as.character(data[[1]]),
      stringsAsFactors = FALSE
    )
  })
  
  observeEvent(input$load_demo, {
    dataset(make_demo_data())
    upload_summary(NULL)
    prediction_values(NULL)
  })
  
  observeEvent(input$demo_file, {
    req(input$demo_file)
    
    uploaded <- tryCatch(
      read_expression_upload(input$demo_file),
      error = function(e) {
        showNotification(
          paste("Upload failed:", conditionMessage(e)),
          type = "error",
          duration = 10
        )
        NULL
      }
    )
    req(uploaded)
    
    dataset(uploaded$data)
    upload_summary(uploaded)
    prediction_values(NULL)
  })
  
  output$sample_selector <- renderUI({
    lookup <- patient_lookup()
    choices <- stats::setNames(lookup$sample_id, lookup$patient_display_id)
    
    tagList(
      selectInput(
        "selected_sample",
        "Select patient",
        choices = choices,
        selected = lookup$sample_id[1]
      ),
      textInput(
        "patient_search",
        "Search patient name",
        placeholder = "Try Patient_1 or Patient_25"
      )
    )
  })
  
  observeEvent(input$patient_search, {
    query <- trimws(input$patient_search)
    if (!nzchar(query)) return()
    
    lookup <- patient_lookup()
    query_lower <- tolower(query)
    
    exact_match <- which(tolower(lookup$patient_display_id) == query_lower)
    if (length(exact_match) == 0) {
      exact_match <- which(tolower(lookup$sample_id) == query_lower)
    }
    
    partial_match <- which(grepl(query_lower, tolower(lookup$patient_display_id), fixed = TRUE))
    if (length(partial_match) == 0) {
      partial_match <- which(grepl(query_lower, tolower(lookup$sample_id), fixed = TRUE))
    }
    
    match_index <- if (length(exact_match) > 0) exact_match[1] else partial_match[1]
    
    if (is.na(match_index)) {
      showNotification(
        paste("No patient found for:", query),
        type = "warning",
        duration = 4
      )
      return()
    }
    
    updateSelectInput(
      session,
      "selected_sample",
      selected = lookup$sample_id[match_index]
    )
  }, ignoreInit = TRUE)
  
  gene_summary <- reactive({
    req(dataset())
    req(input$selected_sample)
    req(input$plot_gene_count)
    
    data <- dataset()
    
    gene_data <- data[, -1, drop = FALSE]
    gene_data <- as.data.frame(lapply(gene_data, function(x) suppressWarnings(as.numeric(x))))
    
    available_genes <- intersect(final_top_genes, colnames(gene_data))
    if (length(available_genes) <= 1) {
      numeric_genes <- names(gene_data)[vapply(gene_data, function(x) any(is.finite(x), na.rm = TRUE), logical(1))]
      available_genes <- head(numeric_genes, 50)
    }
    
    validate(
      need(length(available_genes) > 1, "Not enough numeric gene columns to create the plot.")
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
  
  multiclass_prediction <- reactive({
    req(prediction_values())
    req(dataset())
    req(input$selected_sample)
    
    risk <- prediction_values()$meta
    if (risk < 0.4) {
      return(list(show = FALSE))
    }
    
    if (!isTRUE(multiclass_model$available)) {
      return(list(
        show = TRUE,
        can_score = FALSE,
        message = multiclass_model$message
      ))
    }
    
    data <- dataset()
    selected_index <- which(data[[1]] == input$selected_sample)[1]
    if (is.na(selected_index)) {
      return(list(
        show = TRUE,
        can_score = FALSE,
        message = "Selected sample could not be found in the uploaded dataset."
      ))
    }
    
    gene_data <- data[, -1, drop = FALSE]
    gene_data <- as.data.frame(lapply(gene_data, function(x) suppressWarnings(as.numeric(x))))
    
    panel_genes <- multiclass_model$panel_genes
    present_genes <- intersect(panel_genes, colnames(gene_data))
    missing_genes <- setdiff(panel_genes, colnames(gene_data))
    min_required <- ceiling(length(panel_genes) * 0.80)
    
    if (length(present_genes) < min_required) {
      return(list(
        show = TRUE,
        can_score = FALSE,
        message = paste0(
          "Subtype calculator needs at least ", min_required, " of ",
          length(panel_genes), " panel genes. This sample has ",
          length(present_genes), ". Upload a dataset with real gene symbols from the subtype panel."
        ),
        genes_used = length(present_genes),
        genes_required = length(panel_genes),
        missing_genes = missing_genes
      ))
    }
    
    selected_values <- as.numeric(gene_data[selected_index, present_genes])
    names(selected_values) <- present_genes
    
    grand_mean <- multiclass_model$batch$grand_mean[present_genes]
    sample_mean <- mean(selected_values, na.rm = TRUE)
    adjusted_values <- selected_values - sample_mean + mean(grand_mean, na.rm = TRUE)
    
    x_i <- as.data.frame(matrix(NA_real_, nrow = 1, ncol = length(panel_genes)))
    names(x_i) <- panel_genes
    x_i[1, present_genes] <- adjusted_values
    
    for (gene in panel_genes) {
      if (!is.finite(x_i[[gene]]) || is.na(x_i[[gene]])) {
        x_i[[gene]] <- multiclass_model$imputation$med[[gene]]
      }
    }
    
    x_s <- scale_with_training_params(x_i, multiclass_model$scale)
    subtype_levels <- multiclass_model$fit_rf[["classes"]]
    
    model_probabilities <- predict_multiclass_model_probabilities(
      x_i = x_i,
      x_s = x_s,
      multiclass_model = multiclass_model,
      subtype_levels = subtype_levels
    )
    
    if (length(model_probabilities) == 0) {
      return(list(
        show = TRUE,
        can_score = FALSE,
        message = "Subtype models could not score this sample after preprocessing."
      ))
    }
    
    uncertainty <- calculate_multiclass_uncertainty(
      model_probabilities,
      model_disagreement_reference = multiclass_model$model_disagreement
    )
    probabilities <- sort(uncertainty$average_probabilities, decreasing = TRUE)
    predicted_subtype <- names(probabilities)[1]
    
    importance <- multiclass_model$fit_rf[["importance"]][, 1]
    importance <- importance[panel_genes]
    importance[is.na(importance)] <- 0
    importance_scaled <- importance / max(importance, na.rm = TRUE)
    
    patient_values <- as.numeric(x_i[1, panel_genes])
    names(patient_values) <- panel_genes
    median_values <- multiclass_model$imputation$med[panel_genes]
    change_log2 <- patient_values - median_values
    contribution <- abs(change_log2) * importance_scaled
    
    top_genes <- names(sort(contribution, decreasing = TRUE, na.last = NA))[1:5]
    top_genes <- top_genes[!is.na(top_genes)]
    top_gene_table <- data.frame(
      Gene = top_genes,
      `Change (fold)` = ifelse(
        change_log2[top_genes] >= 0,
        paste0("+", round(2 ^ change_log2[top_genes], 2), "x"),
        paste0("-", round(2 ^ abs(change_log2[top_genes]), 2), "x")
      ),
      Direction = ifelse(change_log2[top_genes] >= 0, "Up", "Down"),
      `Contribution Score` = round(contribution[top_genes], 2),
      Impact = ifelse(contribution[top_genes] >= stats::quantile(contribution, 0.80, na.rm = TRUE), "High", "Moderate"),
      check.names = FALSE
    )
    
    list(
      show = TRUE,
      can_score = TRUE,
      predicted_subtype = predicted_subtype,
      confidence = probabilities[[1]],
      probabilities = probabilities,
      uncertainty = uncertainty,
      top_gene_table = top_gene_table,
      genes_used = length(present_genes),
      genes_required = length(panel_genes),
      missing_gene_count = length(missing_genes),
      uncertainty_rds_name = multiclass_model$uncertainty_rds_name
    )
  })
  
  output$gene_expression_plot_ui <- renderUI({
    gene_count <- as.numeric(input$plot_gene_count)
    plot_height <- max(430, 240 + gene_count * 26)
    plotOutput("gene_expression_plot", height = paste0(plot_height, "px"))
  })
  
  output$multiclass_section <- renderUI({
    pred <- multiclass_prediction()
    
    if (!isTRUE(pred$show)) return(NULL)
    
    if (!isTRUE(pred$can_score)) {
      return(
        div(
          class = "card",
          h3("Molecular Subtype Assessment"),
          div(class = "sample-subtitle", "Displayed because the main risk category is Medium or High."),
          div(class = "warning-box", pred$message)
        )
      )
    }
    
    probability_rows <- lapply(names(pred$probabilities), function(subtype) {
      pct <- pred$probabilities[[subtype]] * 100
      div(
        class = "subtype-prob-row",
        span(subtype),
        div(
          class = "subtype-bar-track",
          div(class = "subtype-bar-fill", style = paste0("width:", round(pct, 1), "%;"))
        ),
        span(paste0(round(pct, 1), "%"))
      )
    })
    
    uncertainty_source <- if (nzchar(pred$uncertainty_rds_name)) {
      paste("Uncertainty calibrated using:", pred$uncertainty_rds_name)
    } else {
      "Uncertainty is estimated from disagreement between saved subtype models."
    }
    
    div(
      class = "card",
      h3("Molecular Subtype Assessment"),
      div(
        class = "sample-subtitle",
        "This second-stage calculator appears only for Medium or High main-risk results."
      ),
      div(
        class = "subtype-result-grid",
        div(
          class = "subtype-primary",
          div(class = "subtype-label", "Most likely subclass"),
          div(class = "subtype-name", pred$predicted_subtype),
          div(class = "subtype-label", "Probability"),
          div(class = "subtype-prob", paste0(round(pred$confidence * 100, 1), "%")),
          div(class = "subtype-pill", paste("Uncertainty:", pred$uncertainty$level))
        ),
        div(
          h4("Subtype Probabilities"),
          probability_rows,
          div(
            class = "subtype-muted",
            paste0(
              "Model agreement based on ", pred$uncertainty$model_count,
              " subtype models. ", uncertainty_source
            )
          )
        )
      ),
      hr(),
      fluidRow(
        column(
          5,
          h4("Gene Coverage / Uncertainty"),
          tableOutput("multiclass_quality_table")
        ),
        column(
          7,
          h4("Top Contributing Genes"),
          DTOutput("multiclass_gene_table")
        )
      )
    )
  })
  
  output$multiclass_quality_table <- renderTable({
    pred <- multiclass_prediction()
    req(isTRUE(pred$can_score))
    
    data.frame(
      Item = c(
        "Panel genes found",
        "Panel genes missing",
        "Prediction confidence",
        "Model disagreement",
        "Low disagreement cutoff",
        "High disagreement cutoff",
        "Uncertainty level"
      ),
      Value = c(
        paste0(pred$genes_used, " / ", pred$genes_required),
        pred$missing_gene_count,
        paste0(round(pred$confidence * 100, 1), "%"),
        round(pred$uncertainty$disagreement, 3),
        round(pred$uncertainty$low_disagreement_cutoff, 3),
        round(pred$uncertainty$high_disagreement_cutoff, 3),
        pred$uncertainty$level
      ),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = FALSE, spacing = "s")
  
  output$multiclass_gene_table <- renderDT({
    pred <- multiclass_prediction()
    req(isTRUE(pred$can_score))
    
    datatable(
      pred$top_gene_table,
      rownames = FALSE,
      options = list(
        dom = "t",
        ordering = FALSE,
        paging = FALSE,
        autoWidth = TRUE
      )
    )
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
    
    upload_text <- if (is.null(upload_summary())) {
      "Current dataset: built-in demo data"
    } else {
      paste(
        "Current dataset:", upload_summary()$source,
        "\nSamples:", upload_summary()$n_samples,
        "\nGenes/features:", upload_summary()$n_genes,
        "\nPlatform:", upload_summary()$platform
      )
    }
    
    paste(
      upload_text,
      "\n\nMain risk model gene check",
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
