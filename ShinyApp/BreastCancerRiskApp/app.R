library(shiny)
library(rsconnect)
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
    file.path("app_data", "models", "multiclass"),
    file.path("Pre-processing", "models", "multiclass"),
    file.path("..", "Pre-processing", "models", "multiclass"),
    file.path("..", "..", "Pre-processing", "models", "multiclass")
  )
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) == 0) return(NULL)
  existing[1]
}

find_risk_model_dir <- function() {
  candidates <- c(
    file.path("app_data", "models"),
    file.path("Pre-processing", "models"),
    file.path("..", "Pre-processing", "models"),
    file.path("..", "..", "Pre-processing", "models")
  )
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) == 0) return(NULL)
  existing[1]
}

find_risk_results_dir <- function() {
  candidates <- c(
    file.path("app_data", "results"),
    file.path("Pre-processing", "results"),
    file.path("..", "Pre-processing", "results"),
    file.path("..", "..", "Pre-processing", "results"),
    file.path("Pre-processing", "result"),
    file.path("..", "Pre-processing", "result"),
    file.path("..", "..", "Pre-processing", "result")
  )
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) == 0) return(NULL)
  existing[1]
}

find_normal_reference_file <- function() {
  candidates <- c(
    file.path("app_data", "normal_reference", "normal_reference.rds")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) return(NULL)
  existing[1]
}

find_normalised_data_file <- function() {
  candidates <- c(
    file.path("Pre-processing", "data", "normalised", "merged_normalised.rds"),
    file.path("..", "Pre-processing", "data", "normalised", "merged_normalised.rds"),
    file.path("..", "..", "Pre-processing", "data", "normalised", "merged_normalised.rds")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) return(NULL)
  existing[1]
}

load_normal_reference <- function(reference_genes, excluded_datasets = "GSE29044") {
  reference_file <- find_normal_reference_file()
  if (!is.null(reference_file)) {
    ref <- readRDS(reference_file)
    if (is.list(ref) && all(c("expr", "genes", "sample_ids", "median", "sd") %in% names(ref))) {
      available_genes <- intersect(reference_genes, ref$genes)
      if (length(available_genes) >= 2) {
        return(list(
          available = TRUE,
          data_file = normalizePath(reference_file),
          expr = ref$expr[available_genes, , drop = FALSE],
          genes = available_genes,
          sample_ids = ref$sample_ids,
          median = ref$median[available_genes],
          sd = ref$sd[available_genes],
          n_samples = length(ref$sample_ids),
          excluded_datasets = if (!is.null(ref$excluded_datasets)) ref$excluded_datasets else excluded_datasets
        ))
      }
    }
  }
  
  data_file <- find_normalised_data_file()
  if (is.null(data_file)) {
    return(list(
      available = FALSE,
      message = "Could not find Pre-processing/data/normalised/merged_normalised.rds."
    ))
  }
  
  dat <- readRDS(data_file)
  if (!all(c("expr", "meta") %in% names(dat))) {
    return(list(
      available = FALSE,
      message = "merged_normalised.rds must contain expr and meta objects."
    ))
  }
  
  normal_idx <- which(
    as.character(dat$meta$Status) == "Normal" &
      !(as.character(dat$meta$dataset) %in% excluded_datasets)
  )
  if (length(normal_idx) == 0) {
    return(list(
      available = FALSE,
      message = "No eligible normal samples found in merged_normalised.rds after excluding external demo datasets."
    ))
  }
  
  normal_sample_ids <- rownames(dat$meta)[normal_idx]
  normal_sample_ids <- intersect(normal_sample_ids, colnames(dat$expr))
  available_genes <- intersect(reference_genes, rownames(dat$expr))
  
  if (length(available_genes) < 2 || length(normal_sample_ids) < 2) {
    return(list(
      available = FALSE,
      message = "Normal reference does not contain enough diagnostic genes or normal samples."
    ))
  }
  
  expr_normal <- dat$expr[available_genes, normal_sample_ids, drop = FALSE]
  normal_median <- apply(expr_normal, 1, median, na.rm = TRUE)
  normal_sd <- apply(expr_normal, 1, stats::sd, na.rm = TRUE)
  normal_sd[!is.finite(normal_sd) | normal_sd == 0] <- 1
  
  list(
    available = TRUE,
    data_file = normalizePath(data_file),
    expr = expr_normal,
    genes = available_genes,
    sample_ids = normal_sample_ids,
    median = normal_median,
    sd = normal_sd,
    n_samples = length(normal_sample_ids),
    excluded_datasets = excluded_datasets
  )
}

load_risk_model <- function() {
  model_dir <- find_risk_model_dir()
  if (is.null(model_dir)) {
    return(list(available = FALSE, message = "Could not find Pre-processing/models."))
  }
  
  model_file <- file.path(model_dir, "elasticnet_full.rds")
  gene_file <- file.path(model_dir, "final_top_genes.rds")
  if (!file.exists(model_file) || !file.exists(gene_file)) {
    return(list(
      available = FALSE,
      message = "Missing elasticnet_full.rds or final_top_genes.rds."
    ))
  }
  fit <- readRDS(model_file)
  genes <- readRDS(gene_file)
  coef_mat <- tryCatch(as.matrix(stats::coef(fit, s = "lambda.min")), error = function(e) NULL)
  coefficients <- stats::setNames(rep(0, length(genes)), genes)
  intercept <- NA_real_
  if (!is.null(coef_mat)) {
    common_coef_genes <- intersect(genes, rownames(coef_mat))
    coefficients[common_coef_genes] <- as.numeric(coef_mat[common_coef_genes, 1])
    if ("(Intercept)" %in% rownames(coef_mat)) {
      intercept <- as.numeric(coef_mat["(Intercept)", 1])
    }
  }
  nonzero_genes <- names(coefficients)[abs(coefficients) > 1e-10]
  
  bootstrap_files <- list.files(
    model_dir,
    pattern = "^elasticnet\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  bootstrap_files <- bootstrap_files[grepl("bootstrap_", bootstrap_files)]
  
  results_dir <- find_risk_results_dir()
  metric_file <- if (!is.null(results_dir)) {
    list.files(results_dir, pattern = "elastic_metric_summary.*\\.rds$", full.names = TRUE)[1]
  } else {
    NA_character_
  }
  
  list(
    available = TRUE,
    model_dir = normalizePath(model_dir),
    fit = fit,
    genes = genes,
    coefficients = coefficients,
    nonzero_genes = nonzero_genes,
    intercept = intercept,
    bootstrap_fits = lapply(bootstrap_files, readRDS),
    metric_summary = if (!is.na(metric_file) && file.exists(metric_file)) readRDS(metric_file) else NULL,
    metric_summary_name = if (!is.na(metric_file) && file.exists(metric_file)) basename(metric_file) else NULL
  )
}

load_multiclass_model <- function() {
  model_dir <- find_multiclass_model_dir()
  if (is.null(model_dir)) {
    return(list(available = FALSE, message = "Could not find Pre-processing/models/multiclass."))
  }
  
  required_files <- c(
    "fit_elasticnet.rds", "panel_genes.rds", "imputation_params.rds",
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
    fit_elasticnet = readRDS(file.path(model_dir, "fit_elasticnet.rds")),
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
    file.path("app_data", "geo-cache"),
    file.path("Pre-processing", "data", "geo-cache"),
    file.path("..", "Pre-processing", "data", "geo-cache"),
    file.path("..", "..", "Pre-processing", "data", "geo-cache")
  )
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) == 0) return(NULL)
  existing[1]
}

find_geo_annotation_file <- function(platform_id) {
  safe_platform <- gsub("[^A-Za-z0-9_.-]", "_", platform_id)
  candidates <- c(
    file.path("app_data", "geo_annotation", paste0(safe_platform, "_probe_map.rds"))
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) return(NULL)
  existing[1]
}

get_upload_gene_panel <- function() {
  genes <- character(0)
  risk_dir <- find_risk_model_dir()
  if (!is.null(risk_dir)) {
    risk_gene_file <- file.path(risk_dir, "final_top_genes.rds")
    if (file.exists(risk_gene_file)) genes <- c(genes, readRDS(risk_gene_file))
  }
  
  multiclass_dir <- find_multiclass_model_dir()
  if (!is.null(multiclass_dir)) {
    multiclass_gene_file <- file.path(multiclass_dir, "panel_genes.rds")
    if (file.exists(multiclass_gene_file)) genes <- c(genes, readRDS(multiclass_gene_file))
  }
  
  unique(as.character(genes))
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

parse_series_matrix_header <- function(datapath) {
  con <- gzfile(datapath, open = "rt")
  on.exit(close(con), add = TRUE)
  
  platform_id <- NA_character_
  series_platform_ids <- character(0)
  sample_platform_ids <- character(0)
  sample_ids <- character(0)
  table_begin_line <- NA_integer_
  line_number <- 0L
  
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    line_number <- line_number + 1L
    
    if (startsWith(line, "!Series_platform_id")) {
      fields <- scan(text = line, what = character(), sep = "\t", quiet = TRUE)
      series_platform_ids <- c(series_platform_ids, gsub('^"|"$', "", fields[-1]))
    }
    
    if (startsWith(line, "!Sample_geo_accession")) {
      fields <- scan(text = line, what = character(), sep = "\t", quiet = TRUE)
      sample_ids <- gsub('^"|"$', "", fields[-1])
    }
    
    if (startsWith(line, "!Sample_platform_id")) {
      fields <- scan(text = line, what = character(), sep = "\t", quiet = TRUE)
      sample_platform_ids <- gsub('^"|"$', "", fields[-1])
    }
    
    if (startsWith(line, "!series_matrix_table_begin")) {
      table_begin_line <- line_number
      break
    }
  }
  
  sample_platform_ids <- unique(sample_platform_ids[nzchar(sample_platform_ids)])
  series_platform_ids <- unique(series_platform_ids[nzchar(series_platform_ids)])
  if (length(sample_platform_ids) == 1) {
    platform_id <- sample_platform_ids[1]
  } else if (length(series_platform_ids) >= 1) {
    platform_id <- series_platform_ids[1]
  }
  
  list(
    platform_id = platform_id,
    sample_ids = sample_ids,
    table_begin_line = table_begin_line
  )
}

read_geo_series_matrix_lightweight <- function(file_info) {
  header <- parse_series_matrix_header(file_info$datapath)
  if (is.na(header$table_begin_line)) {
    stop("Could not find the GEO series matrix expression table.")
  }
  if (is.na(header$platform_id) || !nzchar(header$platform_id)) {
    stop("Could not identify the GEO platform ID.")
  }
  
  annotation_file <- find_geo_annotation_file(header$platform_id)
  if (is.null(annotation_file)) {
    stop(paste0("No lightweight probe annotation is bundled for ", header$platform_id, "."))
  }
  
  probe_map <- readRDS(annotation_file)
  required_genes <- get_upload_gene_panel()
  if (length(required_genes) > 0) {
    probe_map <- rbind(
      probe_map,
      data.frame(
        probe_id = required_genes,
        gene_symbol = required_genes,
        stringsAsFactors = FALSE
      )
    )
    probe_map <- probe_map[probe_map$gene_symbol %in% required_genes, , drop = FALSE]
  }
  if (nrow(probe_map) == 0) {
    stop("No uploaded probes map to the model gene panel.")
  }
  
  probe_to_gene <- split(probe_map$gene_symbol, probe_map$probe_id)
  kept_values <- list()
  kept_genes <- character(0)
  sample_ids <- header$sample_ids
  
  con <- gzfile(file_info$datapath, open = "rt")
  on.exit(close(con), add = TRUE)
  for (i in seq_len(header$table_begin_line)) readLines(con, n = 1, warn = FALSE)
  
  column_line <- readLines(con, n = 1, warn = FALSE)
  if (length(column_line) == 0) stop("The GEO expression table is empty.")
  column_names <- scan(text = column_line, what = character(), sep = "\t", quiet = TRUE)
  column_names <- gsub('^"|"$', "", column_names)
  if (length(sample_ids) == 0) sample_ids <- column_names[-1]
  
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0 || startsWith(line, "!series_matrix_table_end")) break
    
    fields <- scan(text = line, what = character(), sep = "\t", quiet = TRUE)
    if (length(fields) < 2) next
    probe_id <- gsub('^"|"$', "", fields[1])
    gene_symbols <- unname(probe_to_gene[[probe_id]])
    if (is.null(gene_symbols) || length(gene_symbols) == 0) next
    
    values <- suppressWarnings(as.numeric(gsub('^"|"$', "", fields[-1])))
    if (length(values) != length(sample_ids)) next
    
    for (gene_symbol in unique(gene_symbols)) {
      if (is.na(gene_symbol) || !nzchar(gene_symbol)) next
      kept_values[[length(kept_values) + 1]] <- values
      kept_genes <- c(kept_genes, gene_symbol)
    }
  }
  
  if (length(kept_values) == 0) {
    stop("The uploaded GEO file did not contain probes for the model gene panel.")
  }
  
  expr <- do.call(rbind, kept_values)
  rownames(expr) <- kept_genes
  colnames(expr) <- sample_ids
  
  if (stats::quantile(expr, 0.99, na.rm = TRUE) > 100) {
    expr <- log2(expr + 1)
  }
  
  expr_gene <- limma::avereps(expr, ID = rownames(expr))
  expr_table <- as.data.frame(t(expr_gene), check.names = FALSE)
  expr_table <- data.frame(
    sample_id = rownames(expr_table),
    expr_table,
    check.names = FALSE
  )
  
  list(
    data = expr_table,
    source = paste0("GEO series-matrix TXT.GZ lightweight parser (", header$platform_id, ")"),
    n_probes = length(kept_genes),
    n_genes = nrow(expr_gene),
    n_samples = ncol(expr_gene),
    platform = header$platform_id
  )
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
  
  id_candidates <- c(
    "patientid", "patient_id", "patient", "sample_id", "sampleid",
    "sample", "geo_accession", "gsm", "id"
  )
  id_idx <- match(id_candidates, tolower(names(data)), nomatch = 0)
  id_idx <- id_idx[id_idx > 0][1]
  
  if (is.na(id_idx) || is.null(id_idx)) {
    data$sample_id <- paste0("Patient_", sprintf("%03d", seq_len(nrow(data))))
    data <- data[, c(ncol(data), seq_len(ncol(data) - 1)), drop = FALSE]
  } else {
    names(data)[id_idx] <- "sample_id"
    data <- data[, c(id_idx, setdiff(seq_along(data), id_idx)), drop = FALSE]
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
  looks_like_geo <- grepl("\\.txt\\.gz$|series.*matrix", file_name)
  
  if (looks_like_geo) {
    lightweight_error <- NULL
    lightweight_result <- tryCatch(
      read_geo_series_matrix_lightweight(file_info),
      error = function(e) {
        lightweight_error <<- conditionMessage(e)
        warning("Lightweight GEO parser failed: ", lightweight_error)
        NULL
      }
    )
    if (!is.null(lightweight_result)) return(lightweight_result)
    
    stop(
      paste(
        "This GEO series-matrix file could not be processed by the lightweight parser.",
        lightweight_error,
        "Please upload a supported GEO platform file or a preprocessed CSV/TXT table with sample rows and gene columns."
      )
    )
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

format_risk_score <- function(x, digits = 3) {
  if (!is.finite(x)) return("Not available")
  if (x > 0 && x < 10^-digits) {
    return(paste0("<", formatC(10^-digits, format = "f", digits = digits)))
  }
  formatC(round(x, digits), format = "f", digits = digits)
}

align_probability_vector <- function(prob_vec, subtype_levels) {
  out <- setNames(rep(NA_real_, length(subtype_levels)), subtype_levels)
  common <- intersect(names(prob_vec), subtype_levels)
  out[common] <- prob_vec[common]
  out
}

predict_multiclass_model_probabilities <- function(x_i, x_s, multiclass_model, subtype_levels) {
  probs <- list()
  
  elasticnet_prob <- tryCatch({
    pr <- predict(multiclass_model$fit_elasticnet, as.matrix(x_s), type = "response", s = "lambda.min")
    if (length(dim(pr)) == 3) pr <- pr[, , 1]
    if (is.null(dim(pr))) pr <- t(as.matrix(pr))
    align_probability_vector(stats::setNames(as.numeric(pr[1, ]), colnames(pr)), subtype_levels)
  }, error = function(e) NULL)
  if (!is.null(elasticnet_prob)) probs$ElasticNet <- elasticnet_prob
  
  probs
}

calculate_multiclass_uncertainty <- function(model_probabilities, model_disagreement_reference = NULL) {
  prob_stack <- do.call(rbind, model_probabilities)
  avg_probs <- colMeans(prob_stack, na.rm = TRUE)
  avg_probs <- avg_probs / sum(avg_probs, na.rm = TRUE)
  pred_subtype <- names(avg_probs)[which.max(avg_probs)]
  disagreement <- stats::sd(prob_stack[, pred_subtype], na.rm = TRUE)
  if (!is.finite(disagreement) || is.na(disagreement)) disagreement <- 0
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

prepare_risk_features <- function(data, selected_sample, required_genes, important_genes = required_genes, min_fraction = 0.80) {
  selected_index <- which(data[[1]] == selected_sample)[1]
  if (is.na(selected_index)) {
    return(list(can_score = FALSE, message = "Selected sample could not be found in the uploaded dataset."))
  }
  
  gene_data <- data[, -1, drop = FALSE]
  gene_data <- as.data.frame(lapply(gene_data, function(x) suppressWarnings(as.numeric(x))))
  present_genes <- intersect(required_genes, colnames(gene_data))
  missing_genes <- setdiff(required_genes, colnames(gene_data))
  present_important_genes <- intersect(important_genes, colnames(gene_data))
  missing_important_genes <- setdiff(important_genes, colnames(gene_data))
  min_required <- ceiling(length(important_genes) * min_fraction)
  
  if (length(present_important_genes) < min_required) {
    return(list(
      can_score = FALSE,
      message = paste0(
        "Risk calculator needs at least ", min_required, " of ",
        length(important_genes), " non-zero Elastic Net genes. This sample has ",
        length(present_important_genes), "."
      ),
      genes_used = length(present_genes),
      genes_required = length(required_genes),
      missing_genes = missing_genes,
      missing_important_genes = missing_important_genes
    ))
  }
  
  x <- as.data.frame(matrix(NA_real_, nrow = 1, ncol = length(required_genes)))
  names(x) <- required_genes
  x[1, present_genes] <- as.numeric(gene_data[selected_index, present_genes])
  
  for (gene in required_genes) {
    if (is.na(x[[gene]]) || !is.finite(x[[gene]])) {
      if (gene %in% colnames(gene_data)) {
        x[[gene]] <- stats::median(gene_data[[gene]], na.rm = TRUE)
      } else {
        x[[gene]] <- 0
      }
    }
  }
  
  list(
    can_score = TRUE,
    x = x,
    genes_used = length(present_genes),
    genes_required = length(required_genes),
    missing_gene_count = length(missing_genes)
  )
}

predict_elasticnet_risk <- function(x, risk_model) {
  risk <- as.numeric(predict(
    risk_model$fit,
    newx = as.matrix(x),
    s = "lambda.min",
    type = "response"
  ))
  
  bootstrap_scores <- vapply(risk_model$bootstrap_fits, function(fit) {
    as.numeric(predict(
      fit,
      newx = as.matrix(x),
      s = "lambda.min",
      type = "response"
    ))
  }, numeric(1))
  
  if (length(bootstrap_scores) > 0 && any(is.finite(bootstrap_scores))) {
    ci_low <- as.numeric(stats::quantile(bootstrap_scores, 0.10, na.rm = TRUE))
    ci_high <- as.numeric(stats::quantile(bootstrap_scores, 0.90, na.rm = TRUE))
  } else {
    ci_low <- NA_real_
    ci_high <- NA_real_
  }
  
  list(
    risk = risk,
    ci_low = ci_low,
    ci_high = ci_high,
    bootstrap_n = length(bootstrap_scores)
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
      
      hr(),
      
      h4("Patient Selection"),
      uiOutput("sample_selector"),
      
      selectInput(
        "plot_gene_count",
        "Number of genes shown in plot",
        choices = c(
          "10 genes" = 10,
          "20 genes" = 20,
          "50 genes" = 50
        ),
        selected = 10
      ),
      
      actionButton("run_prediction", "Calculate Risk", class = "btn-primary"),
      
      hr(),
      
      h4("Upload Quality Check"),
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
                  div(class = "uncertainty-note", "Estimated prediction range from bootstrap models")
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
                    "This section compares the selected patient with normal reference samples from the training data. GSE29044 is excluded from this reference because it is reserved as the external demo dataset. Boxplots show normal-reference expression; red triangles mark the selected patient.")
              )
          ),
          uiOutput("gene_expression_plot_ui"),
          uiOutput("abnormal_gene_summary"),
          div(class = "gene-value-title", "Selected Patient Gene Expression Summary"),
          DTOutput("gene_value_table")
      ),
      
      div(class = "card",
          h3("Clinical Use Notes"),
          div(class = "warning-box", verbatimTextOutput("warning_box"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  risk_model <- load_risk_model()
  final_top_genes <- if (isTRUE(risk_model$available)) risk_model$genes else paste0("Gene", 1:50)
  normal_reference <- load_normal_reference(final_top_genes)
  display_risk_genes <- if (isTRUE(risk_model$available) && length(risk_model$nonzero_genes) > 0) {
    risk_model$nonzero_genes
  } else {
    final_top_genes
  }
  multiclass_model <- load_multiclass_model()
  
  dataset <- reactiveVal(NULL)
  upload_summary <- reactiveVal(NULL)
  upload_error <- reactiveVal(NULL)
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
  
  observeEvent(input$demo_file, {
    req(input$demo_file)
    
    uploaded <- tryCatch(
      read_expression_upload(input$demo_file),
      error = function(e) {
        upload_error(conditionMessage(e))
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
    upload_error(NULL)
    prediction_values(NULL)
  })
  
  output$sample_selector <- renderUI({
    tagList(
      selectizeInput(
        "selected_sample",
        "Select patient",
        choices = NULL,
        options = list(
          placeholder = "Upload data first, then choose Patient_1...",
          maxOptions = 1000
        )
      ),
      textInput(
        "patient_search",
        "Search patient name",
        placeholder = "Try Patient_1 or Patient_25"
      ),
      div(
        class = "sample-subtitle",
        if (is.null(dataset())) {
          "Upload an expression dataset to activate patient selection."
        } else {
          paste0("Loaded ", nrow(dataset()), " patients/samples from the uploaded file.")
        }
      )
    )
  })
  
  observeEvent(dataset(), {
    data <- dataset()
    if (is.null(data)) {
      updateSelectizeInput(session, "selected_sample", choices = character(0), selected = character(0), server = TRUE)
      return()
    }
    
    lookup <- patient_lookup()
    choices <- stats::setNames(lookup$sample_id, lookup$patient_display_id)
    
    updateSelectizeInput(
      session,
      "selected_sample",
      choices = c("Choose a patient..." = "", choices),
      selected = "",
      server = TRUE
    )
  }, ignoreInit = FALSE)
  
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
    vals <- prediction_values()
    req(vals)
    req(dataset())
    calculated_sample <- vals$sample_id
    validate(need(nzchar(calculated_sample), "Calculate risk to view the selected patient's gene expression profile."))
    req(input$plot_gene_count)
    
    data <- dataset()
    
    gene_data <- data[, -1, drop = FALSE]
    gene_data <- as.data.frame(lapply(gene_data, function(x) suppressWarnings(as.numeric(x))))
    
    validate(
      need(isTRUE(normal_reference$available), normal_reference$message)
    )
    
    available_genes <- Reduce(intersect, list(final_top_genes, colnames(gene_data), normal_reference$genes))
    
    validate(
      need(
        length(available_genes) > 1,
        "Not enough diagnostic panel genes are shared between the uploaded file and the normal reference."
      )
    )
    
    gene_data <- gene_data[, available_genes, drop = FALSE]
    normal_expr <- normal_reference$expr[available_genes, , drop = FALSE]
    
    selected_index <- which(data[[1]] == calculated_sample)[1]
    
    all_patient_values <- as.numeric(gene_data[selected_index, available_genes])
    names(all_patient_values) <- available_genes
    shifted_score <- (all_patient_values - normal_reference$median[available_genes]) / normal_reference$sd[available_genes]
    
    gene_count <- as.numeric(input$plot_gene_count)
    n_high <- ceiling(gene_count / 2)
    n_low <- floor(gene_count / 2)
    
    high_genes <- names(sort(shifted_score, decreasing = TRUE, na.last = NA))
    low_genes <- names(sort(shifted_score, decreasing = FALSE, na.last = NA))
    
    high_genes <- high_genes[seq_len(min(n_high, length(high_genes)))]
    low_genes <- low_genes[seq_len(min(n_low, length(low_genes)))]
    
    # Keep the clinical display order stable:
    # higher-than-reference genes first, then lower-than-reference genes.
    top_genes <- c(high_genes, setdiff(low_genes, high_genes))
    selected_values <- as.numeric(gene_data[selected_index, top_genes])
    selected_coefficients <- if (isTRUE(risk_model$available)) {
      risk_model$coefficients[top_genes]
    } else {
      rep(NA_real_, length(top_genes))
    }
    names(selected_coefficients) <- top_genes
    estimated_contribution <- selected_values * selected_coefficients
    model_direction <- ifelse(
      selected_coefficients > 0,
      "Increases model risk",
      ifelse(selected_coefficients < 0, "Lowers model risk", "No direct model effect")
    )
    
    normal_median <- normal_reference$median[top_genes]
    
    percentile <- sapply(seq_along(top_genes), function(i) {
      gene_values <- as.numeric(normal_expr[top_genes[i], ])
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
      gene_data = as.data.frame(t(normal_reference$expr[available_genes, , drop = FALSE]), check.names = FALSE),
      top_genes = top_genes,
      selected_values = selected_values,
      cohort_median = normal_median,
      percentile = percentile,
      status = status,
      shift_score = shifted_score[top_genes],
      group = ifelse(top_genes %in% high_genes, "Higher than normal reference", "Lower than normal reference"),
      coefficient = selected_coefficients,
      contribution = estimated_contribution,
      model_direction = model_direction,
      n_high = length(high_genes),
      n_low = length(low_genes),
      genes_available = length(available_genes),
      genes_required = length(final_top_genes),
      reference_n = normal_reference$n_samples,
      reference_excluded = paste(normal_reference$excluded_datasets, collapse = ", ")
    )
  })
  
  report_gene_tables <- reactive({
    gs <- gene_summary()
    
    gene_table <- data.frame(
      Gene = gs$top_genes,
      `Patient Expression` = round(gs$selected_values, 4),
      `Normal Reference Median` = round(gs$cohort_median, 4),
      `Percentile vs Normal Reference` = paste0(gs$percentile, "%"),
      `Deviation From Normal Reference` = round(gs$shift_score, 2),
      `Reference Comparison` = gs$group,
      `Model Coefficient` = round(gs$coefficient, 4),
      `Model Direction` = gs$model_direction,
      `Estimated Model Contribution` = round(gs$contribution, 4),
      `Expression Status` = gs$status,
      check.names = FALSE
    )
    
    summary_df <- data.frame(
      Gene = gs$top_genes,
      `Patient Expression` = round(gs$selected_values, 4),
      Percentile = paste0(gs$percentile, "%"),
      Status = gs$status,
      ShiftScore = gs$shift_score,
      Group = gs$group,
      PercentileNumeric = gs$percentile,
      check.names = FALSE
    )
    
    high_genes <- summary_df[summary_df$Group == "Higher than normal reference", ]
    high_genes <- high_genes[order(-high_genes$ShiftScore), ]
    high_genes <- head(high_genes[, c("Gene", "Patient Expression", "Percentile")], 3)
    
    low_genes <- summary_df[summary_df$Group == "Lower than normal reference", ]
    low_genes <- low_genes[order(low_genes$ShiftScore), ]
    low_genes <- head(low_genes[, c("Gene", "Patient Expression", "Percentile")], 3)
    
    list(
      gene_table = gene_table,
      high_genes = high_genes,
      low_genes = low_genes
    )
  })
  
  multiclass_prediction <- reactive({
    vals <- prediction_values()
    req(vals)
    req(dataset())
    calculated_sample <- vals$sample_id
    req(calculated_sample)
    
    risk <- vals$meta
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
    selected_index <- which(data[[1]] == calculated_sample)[1]
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
    subtype_levels <- names(stats::coef(multiclass_model$fit_elasticnet, s = "lambda.min"))
    
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
    
    coef_list <- stats::coef(multiclass_model$fit_elasticnet, s = "lambda.min")
    subtype_coef <- as.matrix(coef_list[[predicted_subtype]])
    importance <- abs(as.numeric(subtype_coef[, 1]))
    names(importance) <- rownames(subtype_coef)
    importance <- importance[setdiff(names(importance), "(Intercept)")]
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
    text(0, -0.27, format_risk_score(risk), cex = 2.5, font = 2, col = "#ff3b30")
    text(
      0,
      -0.47,
      "Range: 0 (Low Risk) – 1 (High Risk)",
      cex = 0.95,
      font = 2,
      col = "#6e6e73"
    )
  })
  
  output$gene_expression_plot <- renderPlot({
    gs <- gene_summary()
    
    ranked_genes <- gs$top_genes
    ranked_values <- gs$selected_values
    
    # Base R horizontal boxplots draw the first variable at the bottom,
    # so reverse the ranking to place the strongest elevation on top.
    plot_genes <- rev(ranked_genes)
    plot_data <- gs$gene_data[, plot_genes, drop = FALSE]
    selected_values <- ranked_values[match(plot_genes, ranked_genes)]
    
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
      main = paste0("Selected Patient Expression Compared With Normal Reference (n = ", gs$reference_n, ")"),
      cex.axis = 0.9,
      cex.lab = 1.1,
      cex.main = 1.2
    )
    
    points(
      x = selected_values,
      y = seq_along(plot_genes),
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
      Status = gs$status,
      ShiftScore = gs$shift_score,
      Group = gs$group
    )
    
    high_genes <- summary_df[summary_df$Group == "Higher than normal reference", ]
    high_genes <- high_genes[order(-high_genes$ShiftScore), ]
    
    low_genes <- summary_df[summary_df$Group == "Lower than normal reference", ]
    low_genes <- low_genes[order(low_genes$ShiftScore), ]
    
    high_items <- if (nrow(high_genes) == 0) {
      list(div(class = "abnormal-item", span("No genes above the normal reference available.")))
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
      list(div(class = "abnormal-item", span("No genes below the normal reference available.")))
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
        div(class = "abnormal-title-high", paste0("Most Elevated vs Normal Reference (Top ", gs$n_high, ")")),
        high_items
      ),
      div(
        class = "abnormal-card",
        div(class = "abnormal-title-low", paste0("Most Reduced vs Normal Reference (Bottom ", gs$n_low, ")")),
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
      `Patient Expression` = round(gs$selected_values, 4),
      `Normal Reference Median` = round(gs$cohort_median, 4),
      `Percentile vs Normal Reference` = paste0(gs$percentile, "%"),
      `Deviation From Normal Reference` = round(gs$shift_score, 2),
      `Reference Comparison` = gs$group,
      `Model Coefficient` = round(gs$coefficient, 4),
      `Model Direction` = gs$model_direction,
      `Estimated Model Contribution` = round(gs$contribution, 4),
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
    dt <- formatStyle(dt, "Patient Expression", fontWeight = "700")
    
    dt
  })
  
  output$gene_check <- renderText({
    if (is.null(dataset())) {
      if (!is.null(upload_error())) {
        return(paste(
          "Upload failed after the file reached the server.",
          upload_error(),
          "Please try a supported GEO TXT.GZ platform or upload a processed CSV/TXT table with patient/sample rows and gene-expression columns.",
          sep = "\n"
        ))
      }
      return(paste(
        "No dataset uploaded.",
        "Please upload a CSV, TXT, or GEO series-matrix TXT.GZ file.",
        "The file must contain patient/sample rows and gene-expression columns for the final Elastic Net model.",
        sep = "\n"
      ))
    }
    
    data <- dataset()
    available_genes <- colnames(data)
    
    found_genes <- intersect(final_top_genes, available_genes)
    missing_genes <- setdiff(final_top_genes, available_genes)
    found_nonzero_genes <- intersect(display_risk_genes, available_genes)
    missing_nonzero_genes <- setdiff(display_risk_genes, available_genes)
    found_reference_genes <- if (isTRUE(normal_reference$available)) {
      intersect(found_genes, normal_reference$genes)
    } else {
      character(0)
    }
    
    upload_text <- if (is.null(upload_summary())) {
      "Current dataset: uploaded expression data"
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
      "Final selected genes:", length(final_top_genes),
      "\nFound selected genes:", length(found_genes),
      "\nMissing selected genes:", length(missing_genes),
      "\nActive model genes required for risk scoring:", length(display_risk_genes),
      "\nFound active model genes:", length(found_nonzero_genes),
      "\nMissing active model genes:", length(missing_nonzero_genes),
      "\nNormal reference samples:", ifelse(isTRUE(normal_reference$available), normal_reference$n_samples, "not available"),
      "\nExcluded from normal reference:", ifelse(isTRUE(normal_reference$available), paste(normal_reference$excluded_datasets, collapse = ", "), "not available"),
      "\nDiagnostic genes available in normal reference and upload:", length(found_reference_genes),
      "\n\nMissing gene list:",
      ifelse(length(missing_nonzero_genes) == 0, "None", paste(missing_nonzero_genes, collapse = ", "))
    )
  })
  
  observeEvent(input$run_prediction, {
    if (is.null(dataset())) {
      showNotification(
        "Please upload an expression dataset before running the calculator.",
        type = "warning",
        duration = 6
      )
      return()
    }
    
    req(input$selected_sample)
    
    if (!isTRUE(risk_model$available)) {
      showNotification(risk_model$message, type = "error", duration = 8)
      return()
    }
    
    prepared <- prepare_risk_features(
      dataset(),
      input$selected_sample,
      risk_model$genes,
      important_genes = risk_model$nonzero_genes
    )
    if (!isTRUE(prepared$can_score)) {
      showNotification(prepared$message, type = "error", duration = 10)
      prediction_values(NULL)
      return()
    }
    
    risk_pred <- predict_elasticnet_risk(prepared$x, risk_model)
    meta_risk <- risk_pred$risk
    ci_low <- if (is.finite(risk_pred$ci_low)) risk_pred$ci_low else NA_real_
    ci_high <- if (is.finite(risk_pred$ci_high)) risk_pred$ci_high else NA_real_
    
    prediction_values(list(
      model = "Elastic Net",
      sample_id = input$selected_sample,
      meta = meta_risk,
      ci_low = ci_low,
      ci_high = ci_high,
      genes_used = prepared$genes_used,
      genes_required = prepared$genes_required,
      missing_gene_count = prepared$missing_gene_count,
      bootstrap_n = risk_pred$bootstrap_n,
      metric_summary_name = risk_model$metric_summary_name
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
    vals <- prediction_values()
    interval_text <- if (is.finite(vals$ci_low) && is.finite(vals$ci_high)) {
      paste0(format_risk_score(vals$ci_low), " - ", format_risk_score(vals$ci_high))
    } else {
      "Not available"
    }
    
    tagList(
      div(class = "uncertainty-pill", interval_text),
      div(
        class = "uncertainty-note",
        paste0("Elastic Net bootstrap models: ", vals$bootstrap_n)
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
      "This dashboard estimates breast cancer likelihood from the uploaded gene-expression profile.",
      "The uploaded file should contain patient/sample IDs and the diagnostic gene panel used by the model.",
      "Gene-expression boxplots compare the selected patient with normal reference samples from the normalized training dataset, excluding GSE29044 because it is reserved as the external demo dataset.",
      "The uncertainty range shows how much the risk estimate varies across saved bootstrap versions of the model.",
      paste0(
        "Internal validation summary loaded from: ",
        ifelse(isTRUE(risk_model$available) && !is.null(risk_model$metric_summary_name),
               risk_model$metric_summary_name,
               "not available")
      ),
      "This tool is intended for research demonstration and clinical decision support only. It should not be used as a standalone diagnostic test without external clinical validation.",
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
      vals <- prediction_values()
      req(vals$sample_id)
      paste0(vals$sample_id, "_risk_report_", Sys.Date(), ".pdf")
    },
    
    content = function(file) {
      vals <- prediction_values()
      req(vals)
      req(vals$sample_id)
      req(gene_summary())
      
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
          patient_id = vals$sample_id,
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
