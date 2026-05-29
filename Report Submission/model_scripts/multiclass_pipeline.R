# Multiclass PAM50 Pipeline — R Script
# Converted from pam50_multiclass_corrected.qmd
# Run this script to generate all RDS outputs for the report.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(GEOquery)
library(limma)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(class)
library(e1071)
library(randomForest)
library(sva)
library(glmnet)
if (requireNamespace("nnet",    quietly = TRUE)) library(nnet)
if (requireNamespace("MASS",    quietly = TRUE)) library(MASS)
select <- dplyr::select
filter <- dplyr::filter
count  <- dplyr::count
if (requireNamespace("xgboost", quietly = TRUE)) library(xgboost)

set.seed(3888)

datasets_multiclass <- c("GSE25066", "GSE45827", "GSE65194", "GSE21653")

subtype_variants <- list(
  "Luminal A"   = c("luminal a","lum a","luma","lum-a","lum_a","luminala"),
  "Luminal B"   = c("luminal b","lum b","lumb","lum-b","lum_b","luminalb"),
  "Basal/TNBC"  = c("basal","basal-like","basal like","basallike","basal_like",
                     "tnbc","triple negative","triple-negative",
                     "triple neg","triplenegative","er-pr-her2-"),
  "HER2"        = c("her2","her2+","her-2","her2-enriched","her2enriched",
                     "her2e","erbb2"),
  "Normal-like" = c("normal-like","normal like","normallike","normal_like",
                    "molecular subtype: normal")
)


# ─────────────────────────────────────────────────────────────────


load_geo_eset <- function(gse_id) {
  geo_cache_dir <- file.path("data", "geo-cache")
  if (!dir.exists(geo_cache_dir)) dir.create(geo_cache_dir, recursive = TRUE)

  gse_matrix_url <- function(id, host = "ftp.ncbi.nlm.nih.gov") {
    num    <- sub("^GSE", "", id)
    prefix <- substr(num, 1, max(1, nchar(num) - 3))
    paste0("https://", host, "/geo/series/GSE", prefix, "nnn/", id,
           "/matrix/", id, "_series_matrix.txt.gz")
  }

  g <- tryCatch(GEOquery::getGEO(gse_id, destdir = geo_cache_dir), error = function(e) e)
  if (!inherits(g, "error")) {
    if (inherits(g, "ExpressionSet")) return(g)
    if (is.list(g) && length(g) >= 1) return(g[[1]])
  }

  destfile <- file.path(geo_cache_dir, paste0(gse_id, "_series_matrix.txt.gz"))
  for (host in c("ftp.ncbi.nlm.nih.gov", "download.ncbi.nlm.nih.gov")) {
    for (attempt in 1:3) {
      try({
        h <- curl::new_handle()
        curl::handle_setopt(h, connecttimeout = 60, timeout = 1800)
        curl::curl_download(gse_matrix_url(gse_id, host), destfile,
                            quiet = FALSE, handle = h)
      }, silent = TRUE)
      if (file.exists(destfile) && file.info(destfile)$size > 0) break
      Sys.sleep(3 * attempt)
    }
  }
  g2 <- GEOquery::getGEO(filename = destfile)
  if (inherits(g2, "ExpressionSet")) return(g2)
  if (is.list(g2) && length(g2) >= 1) return(g2[[1]])
  stop("Could not load ExpressionSet for ", gse_id)
}

needs_log2 <- function(e) {
  v <- as.numeric(e); v <- v[is.finite(v)]
  if (length(v) < 10) return(FALSE)
  qx <- quantile(v, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
  (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
}

safe_log2 <- function(e) {
  pos    <- e[e > 0 & is.finite(e)]
  pseudo <- min(pos, na.rm = TRUE) / 2
  e2 <- e; e2[!is.finite(e2)] <- NA_real_; e2[e2 <= 0] <- pseudo
  log2(e2)
}

quantile_normalise <- function(e) limma::normalizeBetweenArrays(e, method = "quantile")

collapse_probes <- function(eset) {
  e   <- exprs(eset)
  fd  <- fData(eset)
  sym_col <- grep("gene\\s*symbol|^symbol$", colnames(fd),
                  ignore.case = TRUE, value = TRUE)[1]
  gene_sym_raw <- if (!is.na(sym_col)) as.character(fd[, sym_col]) else
    rep(NA_character_, nrow(fd))
  gene_sym <- ifelse(
    is.na(gene_sym_raw) | trimws(gene_sym_raw) == "",
    rownames(e),
    vapply(gene_sym_raw, function(s) {
      s <- trimws(as.character(s))
      if (is.na(s) || s == "") return(NA_character_)
      trimws(strsplit(s, " /// ", fixed = TRUE)[[1]][1])
    }, character(1)))
  gene_sym[is.na(gene_sym) | gene_sym == ""] <-
    rownames(e)[is.na(gene_sym) | gene_sym == ""]
  limma::avereps(e, ID = gene_sym)
}

norm_label_txt <- function(x) {
  x <- tolower(as.character(x)); x[is.na(x)] <- ""
  trimws(gsub("[^a-z0-9]+", " ", x))
}

map_to_subtype <- function(x_raw) {
  x_chr  <- as.character(x_raw)
  x_norm <- norm_label_txt(x_chr)
  if (x_norm == "") return(NA_character_)
  if (grepl("pam50_class", tolower(x_chr), fixed = TRUE)) {
    m <- regexpr("pam50_class:\\s*", x_chr, ignore.case = TRUE)
    if (m > 0) x_norm <- norm_label_txt(
      substring(x_chr, m + attr(m, "match.length")))
  }
  for (canon in names(subtype_variants)) {
    if (any(vapply(subtype_variants[[canon]],
                   function(kw) grepl(kw, x_norm, fixed = TRUE),
                   logical(1))))
      return(canon)
  }
  NA_character_
}

infer_multiclass_label <- function(pMat) {
  cols_try <- c(
    grep("pam50",  colnames(pMat), ignore.case = TRUE, value = TRUE),
    grep("^characteristics_ch", colnames(pMat), value = TRUE),
    setdiff(colnames(pMat), c(
      grep("pam50", colnames(pMat), ignore.case = TRUE, value = TRUE),
      grep("^characteristics_ch", colnames(pMat), value = TRUE)))
  )
  for (cn in cols_try) {
    lab <- vapply(as.character(pMat[[cn]]), map_to_subtype, character(1))
    y   <- droplevels(factor(lab[!is.na(lab)]))
    if (nlevels(y) >= 2) return(list(label = factor(lab), label_col = cn))
  }
  list(label = factor(rep(NA_character_, nrow(pMat))), label_col = NA_character_)
}

normalised_cache_dir <- file.path("../RDS_Files/data")

multiclass_label_overrides <- list(
  GSE21653 = function(p) {
    raw <- tolower(trimws(as.character(p[["characteristics_ch1.12"]])))
    ifelse(grepl("basal",    raw), "Basal/TNBC",
    ifelse(grepl("erbb2",    raw), "HER2",
    ifelse(grepl("luminala", raw), "Luminal A",
    ifelse(grepl("luminalb", raw), "Luminal B",
    ifelse(grepl("normal",   raw), "Normal-like", NA_character_)))))
  }
)

process_dataset_qc <- function(gse_id) {
  rds_path <- file.path(normalised_cache_dir, paste0(gse_id, "_normalised.rds"))

  derive_labels <- function(p) {
    if (gse_id %in% names(multiclass_label_overrides)) {
      message(gse_id, ": using direct label override (characteristics_ch1.12)")
      list(label     = factor(multiclass_label_overrides[[gse_id]](p)),
           label_col = "characteristics_ch1.12")
    } else {
      infer_multiclass_label(p)
    }
  }

  if (file.exists(rds_path)) {
    message(gse_id, ": loading pre-normalised RDS from ", rds_path)
    dat    <- readRDS(rds_path)
    e_gene <- dat$expr
    p      <- dat$meta
    p$dataset_id <- gse_id
    lab    <- derive_labels(p)
    p$label_mc <- as.character(lab$label)
    return(list(id = gse_id, e_raw = e_gene, e_norm = e_gene,
                pMat = p, label_col = lab$label_col))
  }

  message(gse_id, ": no cached RDS found — downloading from GEO")
  eset   <- load_geo_eset(gse_id)
  e      <- exprs(eset)
  if (needs_log2(e)) e <- safe_log2(e)
  e_norm <- quantile_normalise(e)
  exprs(eset) <- e_norm
  e_gene <- collapse_probes(eset)
  p      <- pData(eset)
  p$dataset_id <- gse_id
  lab    <- derive_labels(p)
  p$label_mc <- as.character(lab$label)

  if (!dir.exists(normalised_cache_dir))
    dir.create(normalised_cache_dir, recursive = TRUE)
  saveRDS(list(expr = e_gene, meta = p), rds_path)
  message(gse_id, ": saved normalised RDS to ", rds_path)

  list(id = gse_id, e_raw = e, e_norm = e_gene,
       pMat = p, label_col = lab$label_col)
}


# ─────────────────────────────────────────────────────────────────

load_status <- data.frame(
  dataset = character(), source = character(),
  n_genes = integer(),   n_samples = integer(),
  n_labeled = integer(), subtypes = character(),
  stringsAsFactors = FALSE
)

# ── Load datasets into named list ─────────────────────────────────────────────
qc_pieces <- setNames(
  lapply(datasets_multiclass, function(id) {
    tryCatch(
      process_dataset_qc(id),
      error = function(e) {
        message("FAILED: ", id, " — ", conditionMessage(e))
        NULL
      }
    )
  }),
  datasets_multiclass
)

# ── Remove any failed datasets ────────────────────────────────────────────────
qc_pieces <- list()
for (id in datasets_multiclass) {
  rds_path    <- file.path(normalised_cache_dir, paste0(id, "_normalised.rds"))
  source_used <- if (file.exists(rds_path)) "RDS cache" else "GEO download"
  
  piece <- tryCatch(process_dataset_qc(id), error = function(e) {
    message("FAILED: ", id, " — ", conditionMessage(e)); NULL
  })
  
  if (!is.null(piece)) {
    qc_pieces[[id]] <- piece
    ok <- !is.na(piece$pMat$label_mc) & piece$pMat$label_mc != ""
    load_status <- rbind(load_status, data.frame(
      dataset   = id, source = source_used,
      n_genes   = nrow(piece$e_norm), n_samples = ncol(piece$e_norm),
      n_labeled = sum(ok),
      subtypes  = paste(sort(unique(piece$pMat$label_mc[ok])), collapse = ", "),
      stringsAsFactors = FALSE
    ))
  } else {
    load_status <- rbind(load_status, data.frame(
      dataset = "FAILED", source = source_used,
      n_genes = NA, n_samples = NA, n_labeled = NA, subtypes = "FAILED",
      stringsAsFactors = FALSE
    ))
  }
}

knitr::kable(load_status,
             caption = "Dataset loading summary")

failed <- load_status$source == "FAILED" | is.na(load_status$n_samples)
if (any(failed)) {
  stop("One or more datasets failed to load: ",
       paste(load_status$dataset[failed], collapse = ", "))
}
cat(sprintf("\n✓ All %d datasets loaded successfully.\n", length(qc_pieces)))

# ─────────────────────────────────────────────────────────────────

set.seed(3888)
qc_box_list <- lapply(names(qc_pieces), function(id) {
  p       <- qc_pieces[[id]]
  n_genes <- min(1500L, nrow(p$e_norm))
  gsamp   <- sample(rownames(p$e_norm), n_genes)
  e_sub   <- p$e_norm[gsamp, , drop = FALSE]
  data.frame(
    expr    = as.vector(e_sub),
    sample  = rep(seq_len(ncol(e_sub)), each = nrow(e_sub)),
    dataset = id,
    subtype = rep(p$pMat$label_mc, each = nrow(e_sub)),
    stringsAsFactors = FALSE
  )
})

qc_box_df <- bind_rows(qc_box_list)

sample_med <- qc_box_df %>%
  group_by(dataset, sample) %>%
  summarise(med = median(expr, na.rm = TRUE), .groups = "drop") %>%
  arrange(dataset, med) %>%
  group_by(dataset) %>%
  mutate(sample_ord = row_number()) %>%
  ungroup()

qc_box_df <- qc_box_df %>%
  left_join(sample_med %>% select(dataset, sample, sample_ord),
            by = c("dataset", "sample"))

ggplot(qc_box_df, aes(x = factor(sample_ord), y = expr, colour = subtype)) +
  geom_boxplot(outlier.shape = NA, size = 0.2) +
  facet_wrap(~dataset, scales = "free_x") +
  theme_bw(base_size = 10) +
  labs(x = "Sample (ordered by median)", y = "log2 expression",
       colour = "Subtype",
       title = "Per-sample distributions after within-dataset quantile normalisation") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

qc_box_df %>%
  group_by(dataset) %>%
  summarise(n_samples   = n_distinct(sample),
            median_expr = round(median(expr, na.rm = TRUE), 3),
            IQR_expr    = round(IQR(expr,    na.rm = TRUE), 3),
            q01         = round(quantile(expr, .01, na.rm = TRUE), 3),
            q99         = round(quantile(expr, .99, na.rm = TRUE), 3),
            .groups = "drop") %>%
  knitr::kable(caption = "Per-dataset expression summary (post normalisation)")


# ─────────────────────────────────────────────────────────────────

for (id in names(qc_pieces)) {
  p   <- qc_pieces[[id]]
  ok  <- !is.na(p$pMat$label_mc) & p$pMat$label_mc != ""
  e_s <- p$e_norm[, ok, drop = FALSE]
  y_s <- p$pMat$label_mc[ok]
  e_v <- e_s[apply(e_s, 1, var, na.rm = TRUE) > 0, ]
  if (nrow(e_v) < 5 || ncol(e_v) < 3) { cat(id, ": skipping PCA\n"); next }
  e_v[is.na(e_v)] <- rowMeans(e_v, na.rm = TRUE)[row(e_v)[is.na(e_v)]]
  pca <- prcomp(t(e_v), scale. = TRUE)
  ve  <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
  df  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], subtype = y_s)
  print(ggplot(df, aes(PC1, PC2, colour = subtype)) +
    geom_point(size = 2, alpha = .8) + theme_bw() +
    labs(x = paste0("PC1 (", ve[1], "%)"),
         y = paste0("PC2 (", ve[2], "%)"),
         colour = "Subtype",
         title  = paste0(id, " — PCA after normalisation")))
}


# ─────────────────────────────────────────────────────────────────

all_genes_vis <- Reduce(intersect, lapply(qc_pieces, function(p) rownames(p$e_norm)))
e_merged_vis  <- do.call(cbind, lapply(qc_pieces, function(p)
  p$e_norm[all_genes_vis, , drop = FALSE]))

p_merged_vis <- do.call(rbind, lapply(qc_pieces, function(p) {
  data.frame(dataset = p$id, subtype = p$pMat$label_mc,
             stringsAsFactors = FALSE)
}))

ok_vis    <- !is.na(p_merged_vis$subtype) & p_merged_vis$subtype != ""
e_vis     <- e_merged_vis[, ok_vis, drop = FALSE]
p_vis     <- p_merged_vis[ok_vis, ]
gv        <- apply(e_vis, 1, var, na.rm = TRUE)
top_g     <- names(sort(gv, decreasing = TRUE))[1:min(2000, sum(gv > 0))]
e_vis_sub <- e_vis[top_g, , drop = FALSE]
e_vis_sub[is.na(e_vis_sub)] <- rowMeans(e_vis_sub, na.rm = TRUE)[
  row(e_vis_sub)[is.na(e_vis_sub)]]

make_merged_pca <- function(e_mat, p_df, label) {
  pca <- prcomp(t(e_mat), scale. = TRUE)
  ve  <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
  df  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                    dataset = p_df$dataset, subtype = p_df$subtype)
  p1 <- ggplot(df, aes(PC1, PC2, colour = dataset)) +
    geom_point(size = 1.5, alpha = 0.7) + theme_bw() +
    labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
         colour = "Dataset", title = paste0(label, " — by dataset"))
  p2 <- ggplot(df, aes(PC1, PC2, colour = subtype)) +
    geom_point(size = 1.5, alpha = 0.7) + theme_bw() +
    labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
         colour = "Subtype", title = paste0(label, " — by subtype"))
  list(p1, p2)
}

pl_before <- make_merged_pca(e_vis_sub, p_vis, "Merged — before ComBat")
print(pl_before[[1]]); print(pl_before[[2]])

e_vis_bc <- tryCatch(
  sva::ComBat(dat = e_vis_sub, batch = p_vis$dataset,
              mod = model.matrix(~ factor(p_vis$subtype)),
              par.prior = TRUE, mean.only = FALSE),
  error = function(e) {
    message("ComBat for merged PCA failed: ", conditionMessage(e))
    limma::removeBatchEffect(e_vis_sub, batch = p_vis$dataset,
                             design = model.matrix(~ factor(p_vis$subtype)))
  }
)
pl_after <- make_merged_pca(e_vis_bc, p_vis, "Merged — after ComBat (viz only)")
print(pl_after[[1]]); print(pl_after[[2]])


# ─────────────────────────────────────────────────────────────────

all_genes_mc <- Reduce(intersect, lapply(qc_pieces, function(p) rownames(p$e_norm)))
cat("Shared genes across all multiclass datasets:", length(all_genes_mc), "\n")

e_merged_mc <- do.call(cbind, lapply(qc_pieces, function(p)
  p$e_norm[all_genes_mc, , drop = FALSE]))

na_frac     <- rowMeans(is.na(e_merged_mc))
e_merged_mc <- e_merged_mc[na_frac < 0.3, , drop = FALSE]
cat("Genes after NA filter:", nrow(e_merged_mc), "\n")

p_merged_mc <- do.call(rbind, lapply(qc_pieces, function(p) {
  data.frame(sample_id  = rownames(p$pMat),
             dataset_id = p$id,
             label_mc   = p$pMat$label_mc,
             stringsAsFactors = FALSE)
}))

ok_mc    <- !is.na(p_merged_mc$label_mc) & p_merged_mc$label_mc != ""


y_mc     <- factor(p_merged_mc$label_mc[ok_mc])
e_mc     <- e_merged_mc[, ok_mc, drop = FALSE]
batch_mc <- p_merged_mc$dataset_id[ok_mc]

cat("\nClass counts (merged):\n"); print(table(y_mc))
cat("\nSamples per dataset:\n");   print(table(batch_mc))


# ─────────────────────────────────────────────────────────────────

p_merged_mc %>%
  filter(!is.na(label_mc) & label_mc != "") %>%
  count(dataset_id, label_mc) %>%
  pivot_wider(names_from = label_mc, values_from = n, values_fill = 0) %>%
  knitr::kable(caption = "Subtype counts per dataset")


# ─────────────────────────────────────────────────────────────────


# ── Macro F1 ──────────────────────────────────────────────────────────────────
macro_f1 <- function(truth, pred) {
  truth <- factor(truth); pred <- factor(pred, levels = levels(truth))
  lvls  <- levels(truth)
  f1s <- vapply(lvls, function(k) {
    tp <- sum(pred == k & truth == k, na.rm = TRUE)
    fp <- sum(pred == k & truth != k, na.rm = TRUE)
    fn <- sum(pred != k & truth == k, na.rm = TRUE)
    p  <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
    r  <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
    if (is.na(p) || is.na(r) || (p + r) == 0) NA_real_ else 2 * p * r / (p + r)
  }, numeric(1))
  mean(f1s, na.rm = TRUE)
}

# ── Macro OvR AUC ─────────────────────────────────────────────────────────────
auc_roc_vec <- function(truth01, score) {
  truth01 <- as.integer(truth01); score <- as.numeric(score)
  n_pos <- sum(truth01 == 1); n_neg <- sum(truth01 == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[truth01 == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

macro_auc_ovr <- function(truth, prob_mat) {
  lvls <- colnames(prob_mat)
  aucs <- vapply(lvls, function(k) {
    tk <- as.integer(as.character(truth) == k)
    tk[is.na(tk)] <- 0L
    n_pos <- sum(tk, na.rm = TRUE)
    if (n_pos == 0L || n_pos == length(tk)) return(NA_real_)
    auc_roc_vec(tk, prob_mat[, k])
  }, numeric(1))
  mean(aucs, na.rm = TRUE)
}

# ── Multiclass MCC (Gorodkin) ─────────────────────────────────────────────────
multiclass_mcc <- function(truth, pred) {
  lvls    <- levels(factor(truth))
  truth_f <- factor(as.character(truth), levels = lvls)
  pred_f  <- factor(as.character(pred),  levels = lvls)
  cm      <- table(Actual = truth_f, Predicted = pred_f)
  N       <- sum(cm)
  mcc_num <- N * sum(diag(cm)) - sum(rowSums(cm) * colSums(cm))
  mcc_den <- sqrt((N^2 - sum(colSums(cm)^2)) * (N^2 - sum(rowSums(cm)^2)))
  if (mcc_den == 0) 0 else mcc_num / mcc_den
}

# ── Balanced accuracy ─────────────────────────────────────────────────────────
balanced_accuracy <- function(truth, pred) {
  truth <- factor(truth); pred <- factor(pred, levels = levels(truth))
  lvls  <- levels(truth)
  recalls <- vapply(lvls, function(k) {
    tp <- sum(pred == k & truth == k, na.rm = TRUE)
    fn <- sum(pred != k & truth == k, na.rm = TRUE)
    if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  }, numeric(1))
  mean(recalls, na.rm = TRUE)
}

# ── ComBat: fit on train, apply to test ───────────────────────────────────────
combat_fit_apply <- function(expr_train, batch_train, mod_train,
                              expr_test,  batch_test) {
  if (length(unique(as.character(batch_train))) < 2) {
    message("Single batch in fold — skipping ComBat.")
    return(list(train = expr_train, test = expr_test))
  }
  corrected_train <- tryCatch(
    sva::ComBat(dat = expr_train, batch = batch_train, mod = mod_train,
                par.prior = TRUE, mean.only = FALSE),
    error = function(e) {
      message("ComBat failed, falling back to removeBatchEffect: ",
              conditionMessage(e))
      limma::removeBatchEffect(expr_train, batch = batch_train,
                               design = mod_train)
    }
  )
  batch_train_f <- factor(as.character(batch_train))
  batch_test_f  <- factor(as.character(batch_test),
                           levels = levels(batch_train_f))
  grand_mean    <- rowMeans(expr_train)
  batch_means   <- sapply(levels(batch_train_f), function(b) {
    idx <- which(batch_train_f == b)
    if (length(idx) == 0) return(rep(0, nrow(expr_train)))
    rowMeans(expr_train[, idx, drop = FALSE])
  })
  corrected_test <- expr_test
  for (b in levels(batch_train_f)) {
    idx <- which(as.character(batch_test_f) == b)
    if (length(idx) == 0) next
    corrected_test[, idx] <- expr_test[, idx, drop = FALSE] -
      (batch_means[, b] - grand_mean)
  }
  novel <- which(!(as.character(batch_test) %in% levels(batch_train_f)))
  if (length(novel) > 0) {
    smeans <- colMeans(expr_test[, novel, drop = FALSE])
    corrected_test[, novel] <-
      sweep(expr_test[, novel, drop = FALSE], 2, smeans, "-") + grand_mean
  }
  list(train = corrected_train, test = corrected_test)
}

# ── Imputation & standardisation ─────────────────────────────────────────────
impute_fit   <- function(X) list(med = apply(X, 2, function(v) {
  m <- median(v, na.rm = TRUE); if (is.na(m)) 0 else m }))
impute_apply <- function(X, fit) {
  X2 <- X
  for (j in seq_len(ncol(X2))) {
    idx <- is.na(X2[, j]); if (any(idx)) X2[idx, j] <- fit$med[j]
  }
  X2
}
std_fit      <- function(X) {
  mu <- colMeans(X); sd <- apply(X, 2, sd); sd[sd == 0] <- 1
  list(mu = mu, sd = sd)
}
std_apply    <- function(X, fit) sweep(sweep(X, 2, fit$mu, "-"), 2, fit$sd, "/")

# ── OvR limma gene selection ──────────────────────────────────────────────────
select_genes_ovr <- function(e_mat, y, n_per_class = 50) {
  ok <- !is.na(y) & y != ""; e <- e_mat[, ok, drop = FALSE]
  y  <- factor(y[ok])
  lvls_s <- make.names(levels(y))
  design  <- model.matrix(~0 + y); colnames(design) <- lvls_s
  fit     <- lmFit(e, design)
  genes   <- character(0)
  for (i in seq_along(levels(y))) {
    lvl <- lvls_s[i]; others <- lvls_s[-i]
    cont_str <- paste0(lvl, " - (", paste(others, collapse = "+"),
                       ")/", length(others))
    cont <- makeContrasts(contrasts = cont_str, levels = design)
    fit2 <- eBayes(contrasts.fit(fit, cont))
    genes <- union(genes, rownames(topTable(fit2, n = n_per_class,
                                            sort.by = "P")))
  }
  genes
}

# ── Stratified split ──────────────────────────────────────────────────────────
stratified_split <- function(y, dataset, train_prop = 0.75) {
  dev_idx <- c(); test_idx <- c()
  for (ds in unique(as.character(dataset))) {
    for (cl in unique(as.character(y))) {
      idx <- which(as.character(dataset) == ds & as.character(y) == cl)
      if (length(idx) == 0) next
      idx   <- sample(idx)
      n_dev <- max(1, min(floor(length(idx) * train_prop), length(idx) - 1))
      dev_idx  <- c(dev_idx,  idx[seq_len(n_dev)])
      test_idx <- c(test_idx, idx[(n_dev + 1):length(idx)])
    }
  }
  list(dev = sort(unique(dev_idx)), test = sort(unique(test_idx)))
}


# ─────────────────────────────────────────────────────────────────


run_lodo_multiclass <- function(X, y, batch, n_per_class = 50,
                                 use_weights = TRUE) {
  y     <- factor(y)
  batch <- as.character(batch)
  perf_list  <- list()
  probs_list <- list()

  for (ds in unique(batch)) {
    tryCatch({

    test_id  <- which(batch == ds)
    train_id <- setdiff(seq_along(batch), test_id)
    y_train  <- droplevels(y[train_id])
    y_test   <- factor(as.character(y[test_id]), levels = levels(y_train))
    if (nlevels(y_train) < 2) return(NULL)

    cb <- combat_fit_apply(
      expr_train = t(X[train_id, , drop = FALSE]),
      batch_train = batch[train_id],
      mod_train   = model.matrix(~y_train),
      expr_test   = t(X[test_id, , drop = FALSE]),
      batch_test  = batch[test_id]
    )
    X_tr_full <- t(cb$train); X_te_full <- t(cb$test)

    genes_sel <- select_genes_ovr(t(X_tr_full), y = y_train,
                                  n_per_class = n_per_class)
    sel_g     <- intersect(genes_sel, colnames(X_tr_full))
    if (length(sel_g) < 10) return(NULL)

    X_tr <- X_tr_full[, sel_g]; X_te <- X_te_full[, sel_g]
    imp  <- impute_fit(X_tr)
    X_tr_i <- impute_apply(X_tr, imp); X_te_i <- impute_apply(X_te, imp)
    sc     <- std_fit(X_tr_i)
    X_tr_s <- std_apply(X_tr_i, sc);  X_te_s <- std_apply(X_te_i, sc)

    lvls   <- levels(y_train); n_cl <- length(lvls); n_te <- nrow(X_te_s)
    tab_tr <- table(y_train);  n_tot <- sum(tab_tr)
    cw <- if (use_weights)
      setNames(vapply(lvls, function(k)
        n_tot / (n_cl * max(1, tab_tr[k])), numeric(1)), lvls)
    else setNames(rep(1, n_cl), lvls)
    sw <- cw[as.character(y_train)]

    y_te_chr <- as.character(y_test)

    # ── record(): stores fold metrics + probs ─────────────────────────────────
    # MCC is computed here per fold from predicted class
    record <- function(nm, pred_chr, prob_mat) {
      perf_list[[length(perf_list) + 1]] <<- tibble::tibble(
        fold      = ds,
        model     = nm,
        accuracy  = mean(pred_chr == y_te_chr, na.rm = TRUE),
        macro_f1  = macro_f1(y_te_chr, pred_chr),
        macro_auc = macro_auc_ovr(y_te_chr, prob_mat),
        mcc       = multiclass_mcc(y_te_chr, pred_chr)
      )
      pd <- as.data.frame(prob_mat)
      pd$fold <- ds; pd$model <- nm; pd$truth <- y_te_chr
      probs_list[[length(probs_list) + 1]] <<- pd
    }

    glmnet_probs <- function(fit, X_te, lvls) {
      pr <- predict(fit, as.matrix(X_te), type = "response",
                    s = "lambda.min")[,, 1]
      if (is.null(dim(pr))) pr <- t(as.matrix(pr))
      pr[, intersect(lvls, colnames(pr)), drop = FALSE]
    }

    knn_prob_mat <- function(X_tr, X_te, y_tr, k = 7) {
      pk    <- class::knn(X_tr, X_te, y_tr, k = k, prob = TRUE)
      kw    <- attr(pk, "prob"); pkc <- as.character(pk)
      lvls_k <- levels(y_tr); n_te_k <- nrow(X_te)
      pm <- matrix(0, n_te_k, length(lvls_k),
                   dimnames = list(NULL, lvls_k))
      for (i in seq_len(n_te_k)) {
        pm[i, pkc[i]] <- kw[i]
        oth <- setdiff(lvls_k, pkc[i])
        if (length(oth)) pm[i, oth] <- (1 - kw[i]) / length(oth)
      }
      pm
    }

    sw0 <- rep(1, nrow(X_tr_s))
    cw0 <- setNames(rep(1, n_cl), lvls)

    # ── Unweighted models ─────────────────────────────────────────────────────

    tryCatch({ fit <- e1071::svm(X_tr_s, y_train, kernel = "radial", cost = 1,
                        gamma = 1/ncol(X_tr_s), probability = TRUE)
      ps <- predict(fit, X_te_s, probability = TRUE)
      pm <- attr(ps, "probabilities")[, lvls, drop = FALSE]
      record("SVM_RBF", as.character(ps), pm)
    }, error = function(e) message("SVM_RBF fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- e1071::svm(X_tr_s, y_train, kernel = "linear",
                        cost = 1, probability = TRUE)
      ps <- predict(fit, X_te_s, probability = TRUE)
      pm <- attr(ps, "probabilities")[, lvls, drop = FALSE]
      record("SVM_Linear", as.character(ps), pm)
    }, error = function(e) message("SVM_Linear fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- glmnet::cv.glmnet(as.matrix(X_tr_s), y_train,
                        family = "multinomial", alpha = 1, weights = sw0,
                        standardize = FALSE, type.multinomial = "grouped")
      pm <- glmnet_probs(fit, X_te_s, lvls)
      record("LassoMC", lvls[max.col(pm)], pm)
    }, error = function(e) message("LassoMC fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- glmnet::cv.glmnet(as.matrix(X_tr_s), y_train,
                        family = "multinomial", alpha = 0, weights = sw0,
                        standardize = FALSE, type.multinomial = "grouped")
      pm <- glmnet_probs(fit, X_te_s, lvls)
      record("RidgeMC", lvls[max.col(pm)], pm)
    }, error = function(e) message("RidgeMC fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- glmnet::cv.glmnet(as.matrix(X_tr_s), y_train,
                        family = "multinomial", alpha = 0.5, weights = sw0,
                        standardize = FALSE, type.multinomial = "grouped")
      pm <- glmnet_probs(fit, X_te_s, lvls)
      record("ElasticNetMC", lvls[max.col(pm)], pm)
    }, error = function(e) message("ElasticNetMC fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- nnet::multinom(y_train ~., data = data.frame(y_train, X_tr_s),
                        weights = sw0, trace = FALSE, MaxNWts = 100000)
      pm <- predict(fit, newdata = data.frame(X_te_s), type = "probs")
      if (is.null(dim(pm))) pm <- t(as.matrix(pm))
      pm <- pm[, intersect(lvls, colnames(pm)), drop = FALSE]
      record("Multinom", lvls[max.col(pm)], pm)
    }, error = function(e) message("Multinom fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- e1071::naiveBayes(x = X_tr_s, y = y_train)
      pm <- predict(fit, X_te_s, type = "raw")
      pm <- pm[, intersect(lvls, colnames(pm)), drop = FALSE]
      record("NaiveBayes", colnames(pm)[max.col(pm)], pm)
    }, error = function(e) message("NaiveBayes fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- MASS::lda(x = X_tr_s, grouping = y_train)
      pred_lda <- predict(fit, X_te_s)
      pm <- pred_lda$posterior[, intersect(lvls, colnames(pred_lda$posterior)),
                               drop = FALSE]
      record("LDA", as.character(pred_lda$class), pm)
    }, error = function(e) message("LDA fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ y_int   <- as.integer(y_train) - 1L
      dtrain  <- xgboost::xgb.DMatrix(data = as.matrix(X_tr_s), label = y_int)
      dtest   <- xgboost::xgb.DMatrix(data = as.matrix(X_te_s))
      fit_xgb <- xgboost::xgb.train(
        params = list(objective = "multi:softprob", eval_metric = "mlogloss",
                      num_class = n_cl, eta = 0.05, max_depth = 4,
                      subsample = 0.8, colsample_bytree = 0.8,
                      min_child_weight = 5, nthread = 1),
        data = dtrain, nrounds = 150, verbose = 0)
      prob_raw <- matrix(predict(fit_xgb, dtest), ncol = n_cl, byrow = TRUE)
      colnames(prob_raw) <- lvls
      record("XGBoost", lvls[max.col(prob_raw)], prob_raw)
    }, error = function(e) message("XGBoost fold ", ds, ": ", conditionMessage(e)))

    # ── Weighted models ───────────────────────────────────────────────────────

    tryCatch({ fit <- e1071::svm(X_tr_s, y_train, kernel = "radial", cost = 1,
                        gamma = 1/ncol(X_tr_s), probability = TRUE,
                        class.weights = cw)
      ps <- predict(fit, X_te_s, probability = TRUE)
      pm <- attr(ps, "probabilities")[, lvls, drop = FALSE]
      record("SVM_RBF_w", as.character(ps), pm)
    }, error = function(e) message("SVM_RBF_w fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- e1071::svm(X_tr_s, y_train, kernel = "linear",
                        cost = 1, probability = TRUE, class.weights = cw)
      ps <- predict(fit, X_te_s, probability = TRUE)
      pm <- attr(ps, "probabilities")[, lvls, drop = FALSE]
      record("SVM_Linear_w", as.character(ps), pm)
    }, error = function(e) message("SVM_Linear_w fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- glmnet::cv.glmnet(as.matrix(X_tr_s), y_train,
                        family = "multinomial", alpha = 1, weights = sw,
                        standardize = FALSE, type.multinomial = "grouped")
      pm <- glmnet_probs(fit, X_te_s, lvls)
      record("LassoMC_w", lvls[max.col(pm)], pm)
    }, error = function(e) message("LassoMC_w fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- glmnet::cv.glmnet(as.matrix(X_tr_s), y_train,
                        family = "multinomial", alpha = 0, weights = sw,
                        standardize = FALSE, type.multinomial = "grouped")
      pm <- glmnet_probs(fit, X_te_s, lvls)
      record("RidgeMC_w", lvls[max.col(pm)], pm)
    }, error = function(e) message("RidgeMC_w fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- glmnet::cv.glmnet(as.matrix(X_tr_s), y_train,
                        family = "multinomial", alpha = 0.5, weights = sw,
                        standardize = FALSE, type.multinomial = "grouped")
      pm <- glmnet_probs(fit, X_te_s, lvls)
      record("ElasticNetMC_w", lvls[max.col(pm)], pm)
    }, error = function(e) message("ElasticNetMC_w fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- nnet::multinom(y_train ~., data = data.frame(y_train, X_tr_s),
                        weights = sw, trace = FALSE, MaxNWts = 100000)
      pm <- predict(fit, newdata = data.frame(X_te_s), type = "probs")
      if (is.null(dim(pm))) pm <- t(as.matrix(pm))
      pm <- pm[, intersect(lvls, colnames(pm)), drop = FALSE]
      record("Multinom_w", lvls[max.col(pm)], pm)
    }, error = function(e) message("Multinom_w fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ fit <- e1071::naiveBayes(x = X_tr_s, y = y_train)
      pm <- predict(fit, X_te_s, type = "raw")
      pm <- pm[, intersect(lvls, colnames(pm)), drop = FALSE]
      record("NaiveBayes_w", colnames(pm)[max.col(pm)], pm)
    }, error = function(e) message("NaiveBayes_w fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ prior_w <- cw / sum(cw)
      fit <- MASS::lda(x = X_tr_s, grouping = y_train, prior = prior_w)
      pred_lda <- predict(fit, X_te_s)
      pm <- pred_lda$posterior[, intersect(lvls, colnames(pred_lda$posterior)),
                               drop = FALSE]
      record("LDA_w", as.character(pred_lda$class), pm)
    }, error = function(e) message("LDA_w fold ", ds, ": ", conditionMessage(e)))

    tryCatch({ y_int   <- as.integer(y_train) - 1L
      dtrain  <- xgboost::xgb.DMatrix(data = as.matrix(X_tr_s), label = y_int,
                                       weight = sw)
      dtest   <- xgboost::xgb.DMatrix(data = as.matrix(X_te_s))
      fit_xgb <- xgboost::xgb.train(
        params = list(objective = "multi:softprob", eval_metric = "mlogloss",
                      num_class = n_cl, eta = 0.05, max_depth = 4,
                      subsample = 0.8, colsample_bytree = 0.8,
                      min_child_weight = 5, nthread = 1),
        data = dtrain, nrounds = 150, verbose = 0)
      prob_raw <- matrix(predict(fit_xgb, dtest), ncol = n_cl, byrow = TRUE)
      colnames(prob_raw) <- lvls
      record("XGBoost_w", lvls[max.col(prob_raw)], prob_raw)
    }, error = function(e) message("XGBoost_w fold ", ds, ": ", conditionMessage(e)))
    
    # ── Ensemble: stacked generalisation (ElasticNet + LDA + SVM-RBF) ────────────
    tryCatch({
      n_te   <- nrow(X_te_s)
      n_tr   <- nrow(X_tr_s)
      
      # ── Step 1: inner 3-fold CV on training fold to get OOF base probs ─────────
      # Required to prevent leakage into the meta-learner.
      set.seed(3888)
      inner_folds  <- 3L
      fold_assign  <- sample(rep(seq_len(inner_folds), length.out = n_tr))
      
      get_inner_probs <- function(Xtr, Xte, ytr, lvls) {
        # Locally recomputed weights — avoids NA from missing classes
        ytr      <- droplevels(ytr)
        tab_i    <- table(ytr)
        n_tot_i  <- sum(tab_i)
        n_cl_i   <- nlevels(ytr)
        cw_i     <- setNames(
          vapply(levels(ytr), function(k)
            n_tot_i / (n_cl_i * max(1L, tab_i[k])), numeric(1)),
          levels(ytr))
        sw_i <- cw_i[as.character(ytr)]
        
        out <- list()
        
        # ElasticNet
        tryCatch({
          f  <- glmnet::cv.glmnet(as.matrix(Xtr), ytr,
                                  family = "multinomial", alpha = 0.5,
                                  weights = sw_i, standardize = FALSE,
                                  type.multinomial = "grouped")
          pr <- predict(f, as.matrix(Xte), type = "response",
                        s = "lambda.min")[,,1]
          if (is.null(dim(pr))) pr <- t(as.matrix(pr))
          m  <- pr[, intersect(lvls, colnames(pr)), drop = FALSE]
          rownames(m) <- NULL
          colnames(m) <- paste0("EN_", colnames(m))
          out[["EN"]] <- m
        }, error = function(e) NULL)
        
        # LDA
        tryCatch({
          f  <- MASS::lda(x = Xtr, grouping = ytr)
          pl <- predict(f, Xte)
          m  <- pl$posterior[, intersect(lvls, colnames(pl$posterior)),
                             drop = FALSE]
          rownames(m) <- NULL
          colnames(m) <- paste0("LDA_", colnames(m))
          out[["LDA"]] <- m
        }, error = function(e) NULL)
        
        # SVM RBF
        tryCatch({
          f  <- e1071::svm(Xtr, ytr, kernel = "radial", cost = 1,
                           gamma = 1/ncol(Xtr), probability = TRUE)
          ps <- predict(f, Xte, probability = TRUE)
          m  <- attr(ps, "probabilities")[, intersect(lvls, colnames(
            attr(ps, "probabilities"))), drop = FALSE]
          rownames(m) <- NULL
          colnames(m) <- paste0("SVM_", colnames(m))
          out[["SVM"]] <- m
        }, error = function(e) NULL)
        
        out
      }
      
      # Collect OOF predictions for meta-learner training
      meta_rows <- lapply(seq_len(inner_folds), function(k) {
        i_te  <- which(fold_assign == k)
        i_tr  <- which(fold_assign != k)
        Xtr_k <- X_tr_s[i_tr, , drop = FALSE]
        Xte_k <- X_tr_s[i_te, , drop = FALSE]
        ytr_k <- y_train[i_tr]
        
        probs_k <- get_inner_probs(Xtr_k, Xte_k, ytr_k, lvls)
        if (length(probs_k) == 0) return(NULL)
        
        # Verify all matrices have correct row count
        probs_k <- Filter(function(m) nrow(m) == length(i_te), probs_k)
        if (length(probs_k) == 0) return(NULL)
        
        stacked <- do.call(cbind, probs_k)
        df      <- as.data.frame(stacked, stringsAsFactors = FALSE)
        df$y_inner <- as.character(y_train[i_te])
        df
      })
      
      meta_train_df <- bind_rows(Filter(function(x)
        !is.null(x) && nrow(x) > 0, meta_rows))
      meta_train_df <- meta_train_df[complete.cases(meta_train_df), ]
      
      if (nrow(meta_train_df) < 10)
        stop("meta train too small: ", nrow(meta_train_df), " rows")
      
      y_meta  <- factor(meta_train_df$y_inner, levels = lvls)
      X_meta  <- as.matrix(meta_train_df[,
                                         setdiff(colnames(meta_train_df), "y_inner")])
      
      # ── Step 2: fit meta-learner on OOF stacked probs ──────────────────────────
      tab_meta  <- table(y_meta)
      n_tot_meta <- sum(tab_meta)
      n_cl_meta  <- nlevels(droplevels(y_meta))
      cw_meta   <- setNames(
        vapply(levels(droplevels(y_meta)), function(k)
          n_tot_meta / (n_cl_meta * max(1L, tab_meta[k])), numeric(1)),
        levels(droplevels(y_meta)))
      sw_meta <- cw_meta[as.character(y_meta)]
      
      meta_fit <- glmnet::cv.glmnet(X_meta, y_meta,
                                    family = "multinomial", alpha = 0.5,
                                    weights = sw_meta, standardize = FALSE,
                                    type.multinomial = "grouped")
      
      # ── Step 3: get test-set base model probs (refit on full training fold) ─────
      test_probs <- get_inner_probs(X_tr_s, X_te_s, y_train, lvls)
      
      # Keep only matrices with correct row count
      test_probs <- Filter(function(m) nrow(m) == n_te, test_probs)
      
      if (length(test_probs) == 0)
        stop("no valid base model predictions for test fold")
      
      X_test_meta <- do.call(cbind, test_probs)
      
      # Align columns to what meta-learner was trained on
      shared_cols <- intersect(colnames(X_meta), colnames(X_test_meta))
      if (length(shared_cols) < 2)
        stop("too few shared meta features: ", length(shared_cols))
      
      # Pad any missing columns with 0
      missing_cols <- setdiff(colnames(X_meta), colnames(X_test_meta))
      if (length(missing_cols) > 0) {
        pad <- matrix(0, nrow = n_te, ncol = length(missing_cols),
                      dimnames = list(NULL, missing_cols))
        X_test_meta <- cbind(X_test_meta, pad)
      }
      X_test_meta <- X_test_meta[, colnames(X_meta), drop = FALSE]
      
      # ── Step 4: meta-learner prediction ────────────────────────────────────────
      pr_meta <- predict(meta_fit,
                         as.matrix(X_test_meta),
                         type = "response",
                         s    = "lambda.min")[,,1]
      if (is.null(dim(pr_meta))) pr_meta <- t(as.matrix(pr_meta))
      pm_ens <- pr_meta[, intersect(lvls, colnames(pr_meta)), drop = FALSE]
      
      if (nrow(pm_ens) != n_te)
        stop("ensemble prob matrix row mismatch: ", nrow(pm_ens), " vs ", n_te)
      
      record("Ensemble_Stack", lvls[max.col(pm_ens)], pm_ens)
      message("Ensemble fold ", ds, ": OK")
      
    }, error = function(e) {
      message("Ensemble fold ", ds, " failed: ", conditionMessage(e))
    })

    }, error = function(e) {
      message("LODO fold ", ds, " failed entirely: ", conditionMessage(e))
    })
  }

  if (length(perf_list) == 0)
    warning("run_lodo_multiclass: no folds completed — check data and batch labels")

  list(perf = bind_rows(perf_list), probs = bind_rows(probs_list))
}

X_mc <- t(e_mc)
cat("Running LODO CV (", length(unique(batch_mc)), "folds)...\n")
set.seed(3888)
res_lodo <- run_lodo_multiclass(X_mc, y_mc, batch_mc,
                                n_per_class = 50, use_weights = TRUE)

if (nrow(res_lodo$perf) == 0)
  stop("res_lodo is empty — all LODO folds failed.")

cat("Models evaluated:", paste(unique(res_lodo$perf$model), collapse = ", "), "\n")
cat("Folds completed:", paste(unique(res_lodo$perf$fold), collapse = ", "), "\n")


# ─────────────────────────────────────────────────────────────────


per_class_sens_spec <- function(truth, pred) {
  lvls <- levels(truth)
  bind_rows(lapply(lvls, function(k) {
    tp <- sum(pred == k & truth == k, na.rm = TRUE)
    fn <- sum(pred != k & truth == k, na.rm = TRUE)
    fp <- sum(pred == k & truth != k, na.rm = TRUE)
    tn <- sum(pred != k & truth != k, na.rm = TRUE)
    tibble::tibble(
      subtype     = k,
      sensitivity = if ((tp + fn) == 0) NA_real_ else round(tp/(tp+fn), 3),
      specificity = if ((tn + fp) == 0) NA_real_ else round(tn/(tn+fp), 3),
      n_true      = tp + fn,
      n_pred_pos  = tp + fp
    )
  }))
}

# ── Balanced accuracy per model per fold ──────────────────────────────────────
ba_per_fold <- res_lodo$probs %>%
  group_by(fold, model) %>%
  group_modify(function(df, keys) {
    prob_cols <- setdiff(colnames(df), c("fold", "model", "truth"))
    lv        <- prob_cols[prob_cols %in% levels(factor(df$truth))]
    if (length(lv) < 2) return(tibble::tibble(balanced_acc = NA_real_))
    pred_chr  <- lv[max.col(df[, lv, drop = FALSE])]
    tibble::tibble(balanced_acc = balanced_accuracy(df$truth, pred_chr))
  }) %>%
  ungroup()

ba_summary <- ba_per_fold %>%
  group_by(model) %>%
  summarise(
    ba_mean = mean(balanced_acc, na.rm = TRUE),   # ← mean
    ba_sd   = sd(balanced_acc,   na.rm = TRUE),
    .groups = "drop"
  )

# ── Main scorecard — ALL medians replaced with means ─────────────────────────
scorecard <- res_lodo$perf %>%
  group_by(model) %>%
  summarise(
    macro_auc_mean = mean(macro_auc, na.rm = TRUE),   # ← mean
    macro_auc_sd   = sd(macro_auc,   na.rm = TRUE),
    macro_f1_mean  = mean(macro_f1,  na.rm = TRUE),   # ← mean
    macro_f1_sd    = sd(macro_f1,    na.rm = TRUE),
    accuracy_mean  = mean(accuracy,  na.rm = TRUE),   # ← mean
    accuracy_sd    = sd(accuracy,    na.rm = TRUE),
    mcc_mean       = mean(mcc,       na.rm = TRUE),   # ← MCC added
    mcc_sd         = sd(mcc,         na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(ba_summary, by = "model") %>%
  # Primary sort: MCC mean; secondary: AUC mean
   
  mutate(
    macro_auc    = sprintf("%.3f ± %.3f", macro_auc_mean, macro_auc_sd),
    macro_f1     = sprintf("%.3f ± %.3f", macro_f1_mean,  macro_f1_sd),
    accuracy     = sprintf("%.3f ± %.3f", accuracy_mean,  accuracy_sd),
    mcc          = sprintf("%.3f ± %.3f", mcc_mean,       mcc_sd),
    balanced_acc = sprintf("%.3f ± %.3f", ba_mean,        ba_sd)
  ) %>%
  arrange(desc(balanced_acc),desc(mcc_mean)) %>%
 mutate(
   rank         = row_number()
 ) %>%
  dplyr::select(rank, model, mcc, macro_auc, macro_f1, accuracy, balanced_acc)

knitr::kable(scorecard,
  caption = "Model scorecard — LODO CV. Format: mean ± SD across folds. Sorted by MCC.",
  col.names = c("Rank", "Model", "MCC", "Macro OvR AUC",
                "Macro F1", "Accuracy", "Balanced Accuracy"))

# ── Best model: primary = MCC mean, secondary = AUC mean ─────────────────────
best_model_mc <- res_lodo$perf %>%
  group_by(model) %>%
  summarise(
    m_auc = mean(macro_auc, na.rm = TRUE),
    m_mcc = mean(mcc,       na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    ba_summary %>% select(model, ba_mean),
    by = "model"
  ) %>%
  arrange(desc(ba_mean), desc(m_mcc)) %>%
  slice(1) %>%
  pull(model)

cat("\nBest model (by mean MCC, then mean AUC):", best_model_mc, "\n")

# ── Per-class sens/spec for best model ───────────────────────────────────────
preds_best <- res_lodo$probs %>%
  filter(model == best_model_mc) %>%
  group_modify(function(df, ...) {
    prob_cols <- setdiff(colnames(df), c("fold", "model", "truth"))
    lv        <- prob_cols[prob_cols %in% unique(df$truth)]
    df$pred   <- lv[max.col(df[, lv, drop = FALSE])]
    df
  }) %>%
  ungroup()

all_lvls      <- levels(y_mc)
sens_spec_tbl <- per_class_sens_spec(
  factor(preds_best$truth, levels = all_lvls),
  factor(preds_best$pred,  levels = all_lvls)
)

cat("\nPer-class sensitivity & specificity —", best_model_mc,
    "(pooled LODO):\n")
knitr::kable(sens_spec_tbl,
  caption = paste0("Per-class sensitivity & specificity — ", best_model_mc,
                   " (pooled LODO predictions)"),
  col.names = c("Subtype", "Sensitivity", "Specificity",
                "N true", "N predicted positive"))

# ── Top 50 OvR genes on full batch-corrected dataset ─────────────────────────
cat("\nRunning OvR gene selection on full batch-corrected dataset...\n")

mod_viz  <- model.matrix(~ y_mc)
e_mc_bc  <- tryCatch(
  sva::ComBat(dat = e_mc, batch = batch_mc, mod = mod_viz,
              par.prior = TRUE, mean.only = FALSE),
  error = function(e) limma::removeBatchEffect(e_mc, batch = batch_mc,
                                               design = mod_viz)
)

lvls_mc  <- levels(y_mc)
design_g <- model.matrix(~0 + y_mc)
colnames(design_g) <- make.names(lvls_mc)
fit_g    <- lmFit(e_mc_bc, design_g)

gene_tables <- bind_rows(lapply(lvls_mc, function(k) {
  others   <- make.names(setdiff(lvls_mc, k))
  cont_str <- paste0(make.names(k), " - (",
                     paste(others, collapse = "+"), ")/", length(others))
  cont <- makeContrasts(contrasts = cont_str, levels = design_g)
  fit2 <- eBayes(contrasts.fit(fit_g, cont))
  topTable(fit2, n = 50, sort.by = "P") %>%
    tibble::rownames_to_column("gene") %>%
    mutate(subtype = k, rank = row_number()) %>%
    dplyr::select(subtype, rank, gene, logFC, AveExpr, P.Value, adj.P.Val)
}))

for (k in lvls_mc) {
  knitr::kable(
    gene_tables %>%
      filter(subtype == k) %>%
      mutate(across(c(logFC, AveExpr), round, 3),
             P.Value   = formatC(P.Value,   format = "e", digits = 2),
             adj.P.Val = formatC(adj.P.Val, format = "e", digits = 2)) %>%
      dplyr::select(-subtype),
    caption = paste0("Top 50 OvR genes — ", k)
  ) %>% print()
}

gene_lists  <- setNames(
  lapply(lvls_mc, function(k) gene_tables %>% filter(subtype == k) %>% pull(gene)),
  lvls_mc)
overlap_mat <- outer(lvls_mc, lvls_mc,
  Vectorize(function(a, b) length(intersect(gene_lists[[a]], gene_lists[[b]]))))
rownames(overlap_mat) <- colnames(overlap_mat) <- lvls_mc

cat("\nGene overlap between subtype panels (shared genes in top 50):\n")
knitr::kable(as.data.frame(overlap_mat),
  caption = "Shared genes between subtype top-50 panels.")
cat("\nTotal unique genes:", length(unique(gene_tables$gene)), "\n")


# ─────────────────────────────────────────────────────────────────

library(forcats)

# ── Plot 1: clinician-friendly model comparison (MCC) ─────────────────────────
plot1_df <- res_lodo$perf %>%
  group_by(model) %>%
  summarise(
    auc_mean = mean(macro_auc, na.rm = TRUE),
    auc_sd   = sd(macro_auc,   na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  arrange(desc(auc_mean)) %>%
  mutate(
    rank    = row_number(),
    is_best = model == best_model_mc,
    label   = if_else(is_best, "Selected Classifier",
                      paste0("Alternative ", rank - 1))
  )

ggplot(plot1_df,
       aes(x = auc_mean, y = fct_reorder(label, auc_mean),
           colour = is_best, size = is_best)) +
  geom_linerange(aes(xmin = pmax(auc_mean - auc_sd, 0),
                     xmax = pmin(auc_mean + auc_sd, 1)),
                 linewidth = 0.8, alpha = 0.5) +
  geom_point() +
  geom_text(data = ~filter(.x, is_best),
            aes(label = sprintf("AUC = %.3f", auc_mean)),
            nudge_y = 0.45, fontface = "bold", size = 4.2,
            colour = "#E87722") +
  geom_vline(xintercept = 0.9, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  scale_colour_manual(values = c("TRUE" = "#E87722", "FALSE" = "grey60")) +
  scale_size_manual(values   = c("TRUE" = 5,         "FALSE" = 3)) +
  scale_x_continuous(limits = c(0.45, 1.02),
                     breaks = c(0.5, 0.7, 0.9, 1.0)) +
  theme_bw(base_size = 13) +
  labs(title    = "Classifier Performance Across Independent Datasets",
       subtitle = "Each point = mean AUC; error bars = ±1 SD",
       x = "Discrimination (AUC)  —  higher is better", y = NULL) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

# ── Plot 2: fold-level boxplot for MCC ────────────────────────────────────────
fold_df <- res_lodo$perf %>%
  mutate(
    highlight = if_else(model == best_model_mc, "Best model", "Other"),
    model     = fct_reorder(model, macro_auc, .fun = mean)
  )

ggplot(fold_df, aes(x = model, y = macro_auc, fill = highlight)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, width = 0.6) +
  scale_fill_manual(
    values = c("Best model" = "#E87722", "Other" = "#2166AC"), name = NULL) +
  coord_flip() + theme_bw(base_size = 12) +
  labs(title = "Per-Fold AUC Distribution — All Models",
       x = NULL, y = "Macro OvR AUC (per fold)") +
  theme(legend.position = "top")

# ── Plot 3: held-out sens/spec/BA for best model ──────────────────────────────
mean_sens <- sens_spec_tbl %>%
  filter(!is.na(sensitivity)) %>%
  summarise(m = mean(sensitivity, na.rm = TRUE)) %>% pull(m)
mean_spec <- sens_spec_tbl %>%
  filter(!is.na(specificity)) %>%
  summarise(m = mean(specificity, na.rm = TRUE)) %>% pull(m)
mean_ba   <- ba_summary %>%
  filter(model == best_model_mc) %>% pull(ba_mean)

plot3_df <- tibble::tibble(
  Metric = factor(c("Balanced\nAccuracy", "Sensitivity", "Specificity"),
                  levels = c("Sensitivity", "Specificity", "Balanced\nAccuracy")),
  Score  = c(mean_ba, mean_sens, mean_spec)
)

ggplot(plot3_df, aes(x = Metric, y = Score, fill = Metric)) +
  geom_col(width = 0.5, alpha = 0.9) +
  geom_hline(yintercept = 0.7, linetype = "dashed", colour = "grey40") +
  geom_text(aes(label = round(Score, 2)), vjust = -0.5,
            fontface = "bold", size = 5) +
  scale_fill_manual(values = c("Sensitivity"        = "#2166AC",
                               "Specificity"        = "#35964F",
                               "Balanced\nAccuracy" = "#6A3D9A")) +
  scale_y_continuous(limits = c(0, 1.05)) +
  theme_bw(base_size = 13) +
  labs(title    = paste0(best_model_mc, " — Pooled LODO Performance"),
       subtitle = "Dashed line = 0.70 reference",
       x = "Metric", y = "Score") +
  theme(legend.position = "none")


# ─────────────────────────────────────────────────────────────────

# ── Confusion matrix — best model (pooled LODO predictions) ──────────────────
cm_df <- preds_best %>%
  filter(!is.na(truth), truth != "", truth != "Normal-like") %>%
  mutate(
    truth = factor(truth, levels = setdiff(all_lvls, "Normal-like")),
    pred  = factor(pred,  levels = setdiff(all_lvls, "Normal-like"))
  ) %>%
  count(truth, pred) %>%
  complete(truth, pred, fill = list(n = 0)) %>%
  group_by(truth) %>%
  mutate(
    row_total = sum(n),
    pct       = ifelse(row_total > 0, n / row_total, 0)
  ) %>%
  ungroup()

ggplot(cm_df, aes(x = pred, y = fct_rev(truth), fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(
    label    = ifelse(n > 0, sprintf("%d\n(%.0f%%)", n, pct * 100), "0"),
    colour   = pct > 0.5
  ), size = 3.8, fontface = "bold") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#C6DBEF", "#2166AC"),
    limits  = c(0, 1),
    name    = "Row %"
  ) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"),
                      guide  = "none") +
  scale_x_discrete(position = "top") +
  theme_bw(base_size = 12) +
  labs(
    title    = paste0("Confusion matrix — ", best_model_mc, " (pooled LODO)"),
    x = "Predicted subtype",
    y = "True subtype"
  ) +
  theme(
    axis.text.x      = element_text(angle = 30, hjust = 0),
    panel.grid       = element_blank(),
    legend.position  = "right"
  )


# ─────────────────────────────────────────────────────────────────

res_lodo$perf %>%
  pivot_longer(cols = c(macro_auc, macro_f1, accuracy, mcc),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
    macro_auc = "Macro OvR AUC", macro_f1 = "Macro F1",
    accuracy  = "Accuracy",      mcc      = "MCC")) %>%
  ggplot(aes(x = reorder(model, value, mean), y = value, fill = model)) +
  geom_boxplot(alpha = 0.75, outlier.size = 0.8) +
  facet_wrap(~metric, scales = "free_y") +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = NULL, title = "LODO CV: model comparison") +
  theme(legend.position = "none")


# ─────────────────────────────────────────────────────────────────


split_mc <- stratified_split(y_mc, batch_mc, train_prop = 0.75)
dev_idx  <- split_mc$dev; hold_idx <- split_mc$test

X_dev  <- X_mc[dev_idx,  , drop = FALSE]; y_dev  <- y_mc[dev_idx]
batch_dev  <- batch_mc[dev_idx]
X_hold <- X_mc[hold_idx, , drop = FALSE]; y_hold <- y_mc[hold_idx]
batch_hold <- batch_mc[hold_idx]
cat("Dev:", nrow(X_dev), " Hold-out:", nrow(X_hold), "\n")

res_dev <- run_lodo_multiclass(X_dev, y_dev, batch_dev, n_per_class = 50)

# ── Dev ranking: primary = mean MCC, secondary = mean AUC ─────────────────────
dev_rank <- res_dev$perf %>%
  group_by(model) %>%
  summarise(
    m_mcc = mean(mcc, na.rm = TRUE),
    m_f1  = mean(macro_f1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    ba_per_fold %>%
      group_by(model) %>%
      summarise(ba_mean = mean(balanced_acc, na.rm = TRUE), .groups = "drop"),
    by = "model"
  ) %>%
  arrange(desc(ba_mean), desc(m_mcc)) %>%
  mutate(dev_rank = row_number())

fit_once_mc <- function(X_tr, y_tr, batch_tr, X_te, batch_te,
                         n_per_class = 50) {
  y_tr  <- factor(y_tr); lvls <- levels(y_tr); n_cl <- length(lvls)
  cb    <- combat_fit_apply(t(X_tr), batch_tr, model.matrix(~y_tr),
                            t(X_te), batch_te)
  X_trc <- t(cb$train); X_tec <- t(cb$test)
  gs    <- intersect(select_genes_ovr(t(X_trc), y = y_tr,
                                      n_per_class = n_per_class),
                     colnames(X_trc))
  if (length(gs) < 10) return(NULL)

  X_t2 <- X_trc[, gs]; X_e2 <- X_tec[, gs]
  imp  <- impute_fit(X_t2)
  X_ti <- impute_apply(X_t2, imp); X_ei <- impute_apply(X_e2, imp)
  sc   <- std_fit(X_ti)
  X_ts <- std_apply(X_ti, sc);    X_es <- std_apply(X_ei, sc)

  tab_tr <- table(y_tr); n_tot <- sum(tab_tr); n_te <- nrow(X_es)
  cw  <- setNames(vapply(lvls, function(k)
    n_tot / (n_cl * max(1, tab_tr[k])), numeric(1)), lvls)
  sw  <- cw[as.character(y_tr)]
  sw0 <- rep(1, length(y_tr))

  out <- list()
  gp  <- function(fit) {
    pr <- predict(fit, as.matrix(X_es), type = "response",
                  s = "lambda.min")[,, 1]
    if (is.null(dim(pr))) pr <- t(as.matrix(pr))
    pr[, intersect(lvls, colnames(pr)), drop = FALSE]
  }
  kp <- function() {
    pk  <- class::knn(X_ts, X_es, y_tr, k = 7, prob = TRUE)
    kw  <- attr(pk, "prob"); pkc <- as.character(pk)
    pm  <- matrix(0, n_te, n_cl, dimnames = list(NULL, lvls))
    for (i in seq_len(n_te)) {
      pm[i, pkc[i]] <- kw[i]; oth <- setdiff(lvls, pkc[i])
      if (length(oth)) pm[i, oth] <- (1 - kw[i]) / length(oth)
    }
    pm
  }

  tryCatch({ fit <- e1071::svm(X_ts, y_tr, kernel="radial", cost=1,
    gamma=1/ncol(X_ts), probability=TRUE)
    ps <- predict(fit, X_es, probability=TRUE)
    out[["SVM_RBF"]] <- attr(ps,"probabilities")[,lvls,drop=FALSE]
  }, error=function(e) NULL)
  tryCatch({ fit <- e1071::svm(X_ts, y_tr, kernel="linear", cost=1,
    probability=TRUE)
    ps <- predict(fit, X_es, probability=TRUE)
    out[["SVM_Linear"]] <- attr(ps,"probabilities")[,lvls,drop=FALSE]
  }, error=function(e) NULL)
  tryCatch({ fit <- glmnet::cv.glmnet(as.matrix(X_ts),y_tr,family="multinomial",
    alpha=1,weights=sw0,standardize=FALSE,type.multinomial="grouped")
    out[["LassoMC"]] <- gp(fit) }, error=function(e) NULL)
  tryCatch({ fit <- glmnet::cv.glmnet(as.matrix(X_ts),y_tr,family="multinomial",
    alpha=0,weights=sw0,standardize=FALSE,type.multinomial="grouped")
    out[["RidgeMC"]] <- gp(fit) }, error=function(e) NULL)
  tryCatch({ fit <- glmnet::cv.glmnet(as.matrix(X_ts),y_tr,family="multinomial",
    alpha=0.5,weights=sw0,standardize=FALSE,type.multinomial="grouped")
    out[["ElasticNetMC"]] <- gp(fit) }, error=function(e) NULL)
  tryCatch({ fit <- nnet::multinom(y_tr~.,data=data.frame(y_tr,X_ts),
    weights=sw0,trace=FALSE,MaxNWts=100000)
    pm <- predict(fit,newdata=data.frame(X_es),type="probs")
    if (is.null(dim(pm))) pm <- t(as.matrix(pm))
    out[["Multinom"]] <- pm[,intersect(lvls,colnames(pm)),drop=FALSE]
  }, error=function(e) NULL)
  tryCatch({ fit <- e1071::naiveBayes(x=X_ts,y=y_tr)
    pm <- predict(fit,X_es,type="raw")
    out[["NaiveBayes"]] <- pm[,intersect(lvls,colnames(pm)),drop=FALSE]
  }, error=function(e) NULL)
  tryCatch({ fit <- MASS::lda(x=X_ts,grouping=y_tr)
    pl <- predict(fit,X_es)
    out[["LDA"]] <- pl$posterior[,intersect(lvls,colnames(pl$posterior)),drop=FALSE]
  }, error=function(e) NULL)
  tryCatch({ y_int <- as.integer(y_tr)-1L
    dtr <- xgboost::xgb.DMatrix(data=as.matrix(X_ts),label=y_int)
    dte <- xgboost::xgb.DMatrix(data=as.matrix(X_es))
    fit <- xgboost::xgb.train(params=list(objective="multi:softprob",
      eval_metric="mlogloss",num_class=n_cl,eta=0.05,max_depth=4,
      subsample=0.8,colsample_bytree=0.8,min_child_weight=5,nthread=1),
      data=dtr,nrounds=150,verbose=0)
    pm <- matrix(predict(fit,dte),ncol=n_cl,byrow=TRUE); colnames(pm)<-lvls
    out[["XGBoost"]] <- pm }, error=function(e) NULL)

  tryCatch({ fit <- e1071::svm(X_ts,y_tr,kernel="radial",cost=1,
    gamma=1/ncol(X_ts),probability=TRUE,class.weights=cw)
    ps <- predict(fit,X_es,probability=TRUE)
    out[["SVM_RBF_w"]] <- attr(ps,"probabilities")[,lvls,drop=FALSE]
  }, error=function(e) NULL)
  tryCatch({ fit <- e1071::svm(X_ts,y_tr,kernel="linear",cost=1,
    probability=TRUE,class.weights=cw)
    ps <- predict(fit,X_es,probability=TRUE)
    out[["SVM_Linear_w"]] <- attr(ps,"probabilities")[,lvls,drop=FALSE]
  }, error=function(e) NULL)
  tryCatch({ fit <- glmnet::cv.glmnet(as.matrix(X_ts),y_tr,family="multinomial",
    alpha=1,weights=sw,standardize=FALSE,type.multinomial="grouped")
    out[["LassoMC_w"]] <- gp(fit) }, error=function(e) NULL)
  tryCatch({ fit <- glmnet::cv.glmnet(as.matrix(X_ts),y_tr,family="multinomial",
    alpha=0,weights=sw,standardize=FALSE,type.multinomial="grouped")
    out[["RidgeMC_w"]] <- gp(fit) }, error=function(e) NULL)
  tryCatch({ fit <- glmnet::cv.glmnet(as.matrix(X_ts),y_tr,family="multinomial",
    alpha=0.5,weights=sw,standardize=FALSE,type.multinomial="grouped")
    out[["ElasticNetMC_w"]] <- gp(fit) }, error=function(e) NULL)
  tryCatch({ fit <- nnet::multinom(y_tr~.,data=data.frame(y_tr,X_ts),
    weights=sw,trace=FALSE,MaxNWts=100000)
    pm <- predict(fit,newdata=data.frame(X_es),type="probs")
    if (is.null(dim(pm))) pm <- t(as.matrix(pm))
    out[["Multinom_w"]] <- pm[,intersect(lvls,colnames(pm)),drop=FALSE]
  }, error=function(e) NULL)
  tryCatch({ fit <- e1071::naiveBayes(x=X_ts,y=y_tr)
    pm <- predict(fit,X_es,type="raw")
    out[["NaiveBayes_w"]] <- pm[,intersect(lvls,colnames(pm)),drop=FALSE]
  }, error=function(e) NULL)
  tryCatch({ prior_w <- cw/sum(cw)
    fit <- MASS::lda(x=X_ts,grouping=y_tr,prior=prior_w)
    pl  <- predict(fit,X_es)
    out[["LDA_w"]] <- pl$posterior[,intersect(lvls,colnames(pl$posterior)),drop=FALSE]
  }, error=function(e) NULL)
  tryCatch({ y_int <- as.integer(y_tr)-1L
    dtr <- xgboost::xgb.DMatrix(data=as.matrix(X_ts),label=y_int,weight=sw)
    dte <- xgboost::xgb.DMatrix(data=as.matrix(X_es))
    fit <- xgboost::xgb.train(params=list(objective="multi:softprob",
      eval_metric="mlogloss",num_class=n_cl,eta=0.05,max_depth=4,
      subsample=0.8,colsample_bytree=0.8,min_child_weight=5,nthread=1),
      data=dtr,nrounds=150,verbose=0)
    pm <- matrix(predict(fit,dte),ncol=n_cl,byrow=TRUE); colnames(pm)<-lvls
    out[["XGBoost_w"]] <- pm }, error=function(e) NULL)

  Filter(Negate(is.null), out)
}

hold_probs <- fit_once_mc(X_dev, y_dev, batch_dev, X_hold, batch_hold)

if (!is.null(hold_probs)) {
  lvls_hold <- levels(y_dev)

  single_hold <- bind_rows(lapply(names(hold_probs), function(nm) {
    tryCatch({
      pm  <- hold_probs[[nm]]
      if (is.null(pm) || ncol(pm) == 0) return(NULL)
      lv  <- intersect(lvls_hold, colnames(pm))
      if (length(lv) < 2) return(NULL)
      pred_chr <- lv[max.col(pm[, lv, drop = FALSE])]
      tibble::tibble(
        type      = "Single model", id = nm,
        accuracy  = mean(pred_chr == as.character(y_hold), na.rm = TRUE),
        macro_f1  = macro_f1(as.character(y_hold), pred_chr),
        macro_auc = macro_auc_ovr(as.character(y_hold), pm[, lv, drop = FALSE]),
        mcc       = multiclass_mcc(as.character(y_hold), pred_chr)
      )
    }, error = function(e) NULL)
  }))

  if (nrow(single_hold) == 0)
    stop("No models produced valid hold-out predictions")

  # ── Sort by MCC mean ──────────────────────────────────────────────────────
  single_hold <- single_hold %>% arrange(desc(mcc))

  knitr::kable(
    single_hold %>% mutate(across(c(macro_auc, macro_f1, accuracy, mcc),
                                  round, 3)),
    caption = "Hold-out: single model performance ranked by MCC"
  )

  single_hold %>%
    mutate(across(c(macro_auc, macro_f1, accuracy, mcc), round, 3)) %>%
    pivot_longer(cols = c(macro_auc, macro_f1, accuracy, mcc),
                 names_to = "metric", values_to = "value") %>%
    mutate(metric   = recode(metric, macro_auc = "Macro OvR AUC",
                             macro_f1 = "Macro F1", accuracy = "Accuracy",
                             mcc = "MCC"),
           weighted = grepl("_w$", id)) %>%
    ggplot(aes(x = reorder(id, value), y = value, fill = weighted)) +
    geom_col(alpha = 0.85) +
    geom_text(aes(label = sprintf("%.3f", value)), hjust = -0.1, size = 2.8) +
    facet_wrap(~metric, scales = "free_x") +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = "grey50", "TRUE" = "#D7191C"),
                      labels = c("FALSE" = "Unweighted", "TRUE" = "Weighted")) +
    theme_bw(base_size = 9) +
    labs(x = NULL, y = NULL, fill = NULL,
         title = "Hold-out performance — all models") +
    theme(legend.position = "bottom")
}


# ─────────────────────────────────────────────────────────────────

if (!is.null(hold_probs) && length(hold_probs) >= 2) {
  lvls_u     <- lvls_hold
  sel_mods_u <- dev_rank %>% slice_head(n = min(4, nrow(dev_rank))) %>%
    pull(model) %>% intersect(names(hold_probs))

  prob_stack <- array(NA_real_,
    dim = c(nrow(X_hold), length(lvls_u), length(sel_mods_u)),
    dimnames = list(NULL, lvls_u, sel_mods_u))
  for (nm in sel_mods_u) {
    pm <- hold_probs[[nm]]; lv <- intersect(lvls_u, colnames(pm))
    prob_stack[, lv, nm] <- pm[, lv]
  }

  avg_stack  <- apply(prob_stack, c(1, 2), mean, na.rm = TRUE)
  pred_class <- lvls_u[max.col(avg_stack)]

  pred_class_idx <- match(pred_class, lvls_u)
  disagreement   <- vapply(seq_len(nrow(X_hold)), function(i)
    sd(prob_stack[i, pred_class_idx[i], ], na.rm = TRUE), numeric(1))

  entropy <- apply(avg_stack, 1, function(p) {
    p <- p[p > 0]; -sum(p * log(p))
  })

  uncertainty_df <- tibble::tibble(
    sample_idx   = seq_len(nrow(X_hold)),
    true_subtype = as.character(y_hold),
    pred_subtype = pred_class,
    correct      = pred_class == as.character(y_hold),
    confidence   = apply(avg_stack, 1, max),
    disagreement = disagreement,
    entropy      = entropy
  )

  p1 <- ggplot(uncertainty_df,
       aes(x = confidence, y = disagreement, colour = correct)) +
  geom_jitter(alpha = 0.65, size = 2, width = 0.008, height = 0.003,
              seed = 42) +
  scale_colour_manual(values = c("TRUE" = "#41AB5D", "FALSE" = "#D7191C"),
                      labels = c("TRUE" = "Correct", "FALSE" = "Incorrect")) +
  theme_bw() +
  labs(x = "Ensemble confidence (max P)",
       y = "Model disagreement (SD of P)",
       colour = NULL,
       title    = "Confidence vs model disagreement",
       subtitle = "Low-right = high confidence, models agree; top-left = uncertain") +
  theme(legend.position = "right")

  p2 <- ggplot(uncertainty_df,
               aes(x = true_subtype, y = disagreement, fill = true_subtype)) +
    geom_boxplot(alpha = 0.75, outlier.size = 0.8) +
    theme_bw() +
    labs(x = "True subtype", y = "Model disagreement (SD)",
         title = "Which subtypes are hardest to predict consistently?") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))

  print(p1); print(p2)

  threshold_u <- quantile(uncertainty_df$disagreement, 0.75, na.rm = TRUE)
  knitr::kable(
    uncertainty_df %>%
      filter(disagreement > threshold_u) %>%
      arrange(desc(disagreement)) %>%
      mutate(across(c(confidence, disagreement, entropy), round, 3)) %>%
      head(15),
    caption = "High-uncertainty samples — models disagree most on these"
  )
}


# ─────────────────────────────────────────────────────────────────

set.seed(3888)
B_boot_mc <- 200

boot_model_name <- best_model_mc
boot_alpha <- if (grepl("Lasso",      boot_model_name)) 1   else
              if (grepl("Ridge",      boot_model_name)) 0   else
              if (grepl("ElasticNet", boot_model_name)) 0.5 else NA

cat("Bootstrap will refit:", boot_model_name, "\n")
if (!is.na(boot_alpha)) cat("  glmnet alpha:", boot_alpha, "\n")

run_bootstrap_mc <- function(X_tr, y_tr, batch_tr, X_te, batch_te,
                              B = 200, n_per_class = 50,
                              model_name = "LassoMC_w", alpha = 1) {
  y_tr  <- factor(y_tr); lvls <- levels(y_tr); n_cl <- length(lvls)
  n_te  <- nrow(X_te)

  cb    <- combat_fit_apply(t(X_tr), batch_tr, model.matrix(~y_tr),
                            t(X_te), batch_te)
  X_trc <- t(cb$train); X_tec <- t(cb$test)

  gs <- tryCatch({
    sel <- select_genes_ovr(t(X_trc), y = y_tr, n_per_class = n_per_class)
    intersect(sel, colnames(X_trc))
  }, error = function(e) character(0))
  if (length(gs) < 10) { message("Too few genes for bootstrap"); return(NULL) }

  X_t2 <- X_trc[, gs]; X_e2 <- X_tec[, gs]
  imp  <- impute_fit(X_t2)
  X_ti <- impute_apply(X_t2, imp); X_ei <- impute_apply(X_e2, imp)
  sc   <- std_fit(X_ti)
  X_ts <- std_apply(X_ti, sc);    X_es <- std_apply(X_ei, sc)

  n_tr   <- nrow(X_ts)
  tab_tr <- table(y_tr); n_tot <- sum(tab_tr)
  use_w  <- grepl("_w$", model_name)
  cw     <- setNames(vapply(lvls, function(k)
    n_tot / (n_cl * max(1, tab_tr[k])), numeric(1)), lvls)
  sw_full <- if (use_w) cw[as.character(y_tr)] else rep(1, n_tr)

  boot_arr <- array(NA_real_, dim = c(n_te, n_cl, B),
                    dimnames = list(NULL, lvls, NULL))

  if (!is.na(alpha)) {
    cv_once <- tryCatch(
      glmnet::cv.glmnet(as.matrix(X_ts), y_tr, family = "multinomial",
                        alpha = alpha, weights = sw_full, standardize = FALSE,
                        type.multinomial = "grouped"),
      error = function(e) NULL
    )
    if (is.null(cv_once)) { message("cv.glmnet setup failed"); return(NULL) }
    lambda_fixed_mc <- cv_once$lambda.min
    cat(sprintf("  Fixed lambda: %.5f\n", lambda_fixed_mc))

    for (b in seq_len(B)) {
      idx_b <- sample(n_tr, n_tr, replace = TRUE)
      y_b   <- droplevels(y_tr[idx_b])
      X_b   <- X_ts[idx_b, , drop = FALSE]
      if (nlevels(y_b) < 2) next
      tab_b  <- table(y_b); n_tot_b <- sum(tab_b); lvls_b <- levels(y_b)
      sw_b   <- if (use_w)
        setNames(vapply(lvls_b, function(k)
          n_tot_b / (length(lvls_b) * max(1, tab_b[k])), numeric(1)),
          lvls_b)[as.character(y_b)]
      else rep(1, length(y_b))
      lf <- tryCatch(
        glmnet::glmnet(as.matrix(X_b), y_b, family = "multinomial",
                       alpha = alpha, weights = sw_b, standardize = FALSE,
                       lambda = lambda_fixed_mc, type.multinomial = "grouped"),
        error = function(e) NULL
      )
      if (is.null(lf)) next
      pr <- predict(lf, as.matrix(X_es), type = "response",
                    s = lambda_fixed_mc)[,, 1]
      if (is.null(dim(pr))) pr <- t(as.matrix(pr))
      lv_common <- intersect(lvls, colnames(pr))
      boot_arr[, lv_common, b] <- pr[, lv_common]
    }
  } else if (grepl("Multinom", model_name)) {
    for (b in seq_len(B)) {
      idx_b <- sample(n_tr, n_tr, replace = TRUE)
      y_b   <- droplevels(y_tr[idx_b])
      X_b   <- X_ts[idx_b, , drop = FALSE]
      if (nlevels(y_b) < 2) next
      tab_b  <- table(y_b); n_tot_b <- sum(tab_b); lvls_b <- levels(y_b)
      sw_b   <- if (use_w)
        setNames(vapply(lvls_b, function(k)
          n_tot_b / (length(lvls_b) * max(1, tab_b[k])), numeric(1)),
          lvls_b)[as.character(y_b)]
      else rep(1, length(y_b))
      mn <- tryCatch(
        nnet::multinom(y_b ~., data = data.frame(y_b, X_b),
                       weights = sw_b, trace = FALSE, MaxNWts = 100000),
        error = function(e) NULL
      )
      if (is.null(mn)) next
      pm <- predict(mn, newdata = data.frame(X_es), type = "probs")
      if (is.null(dim(pm))) pm <- t(as.matrix(pm))
      lv_common <- intersect(lvls, colnames(pm))
      boot_arr[, lv_common, b] <- pm[, lv_common]
    }
  } else {
    message("Bootstrap not implemented for model: ", model_name)
    return(NULL)
  }

  boot_arr
}

cat("Running", B_boot_mc, "bootstrap iterations using:", boot_model_name, "\n")
boot_arr_mc <- run_bootstrap_mc(
  X_dev, y_dev, batch_dev, X_hold, batch_hold,
  B = B_boot_mc, n_per_class = 50,
  model_name = boot_model_name, alpha = boot_alpha
)

if (!is.null(boot_arr_mc)) {
  lvls_b    <- dimnames(boot_arr_mc)[[2]]
  prob_mean <- apply(boot_arr_mc, c(1, 2), mean, na.rm = TRUE)
  prob_lo   <- apply(boot_arr_mc, c(1, 2), quantile, 0.025, na.rm = TRUE)
  prob_hi   <- apply(boot_arr_mc, c(1, 2), quantile, 0.975, na.rm = TRUE)
  ci_width  <- prob_hi - prob_lo

  pred_subtype_boot <- lvls_b[max.col(prob_mean)]
  confidence_boot   <- apply(prob_mean, 1, max)
  ci_of_pred        <- vapply(seq_len(nrow(prob_mean)), function(i)
    ci_width[i, pred_subtype_boot[i]], numeric(1))

  boot_risk_mc <- tibble::tibble(
    sample_id    = seq_len(nrow(X_hold)),
    true_subtype = as.character(y_hold),
    pred_subtype = pred_subtype_boot,
    correct      = pred_subtype_boot == as.character(y_hold),
    confidence   = round(confidence_boot, 3),
    ci_width_pred = round(ci_of_pred, 3),
    uncertainty_band = case_when(
      ci_of_pred <= 0.10 ~ "Narrow (confident)",
      ci_of_pred <= 0.20 ~ "Moderate",
      TRUE               ~ "Wide (uncertain — review)"
    )
  )

  for (k in lvls_b) {
    boot_risk_mc[[paste0("P_", k)]]     <- round(prob_mean[, k], 3)
    boot_risk_mc[[paste0("CI_lo_", k)]] <- round(prob_lo[,  k], 3)
    boot_risk_mc[[paste0("CI_hi_", k)]] <- round(prob_hi[,  k], 3)
  }

  cat("\nUncertainty band distribution:\n")
  print(table(boot_risk_mc$uncertainty_band))
  cat("\nAccuracy (bootstrap mean prediction):",
      round(mean(boot_risk_mc$correct, na.rm = TRUE), 3), "\n")

  knitr::kable(
    boot_risk_mc %>%
      arrange(desc(confidence)) %>%
      select(sample_id, true_subtype, pred_subtype, correct,
             confidence, ci_width_pred, uncertainty_band) %>%
      head(20),
    caption = paste0("Per-patient PAM50 predictions with bootstrap 95% CI (B=",
                     B_boot_mc, ").")
  )

  boot_risk_mc %>%
    arrange(desc(confidence)) %>%
    slice_head(n = 40) %>%
    mutate(sample_id = factor(sample_id, levels = rev(sample_id))) %>%
    ggplot(aes(x = confidence, y = sample_id, colour = correct)) +
    geom_errorbarh(aes(xmin = confidence - ci_width_pred / 2,
                       xmax = confidence + ci_width_pred / 2),
                   height = 0.4, alpha = 0.6) +
    geom_point(size = 2) +
    scale_colour_manual(values = c("TRUE" = "#41AB5D", "FALSE" = "#D7191C"),
                        labels = c("TRUE" = "Correct", "FALSE" = "Incorrect")) +
    theme_bw(base_size = 9) +
    labs(x = "Predicted subtype confidence [point = mean, bar = 95% CI]",
         y = "Hold-out patient", colour = NULL,
         title = paste0("Per-patient subtype confidence (B=", B_boot_mc, ")"),
         subtitle = "Wide bars = model uncertain. Dashed = 0.5 threshold.") +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "bottom")

  ggplot(boot_risk_mc,
         aes(x = true_subtype, y = ci_width_pred, fill = true_subtype)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.8) +
    theme_bw() +
    labs(x = "True subtype",
         y = "Bootstrap 95% CI width on predicted subtype",
         title = "Which subtypes are hardest to predict with confidence?") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))

  cat(sprintf("\nMedian CI width: %.3f  |  Mean CI width: %.3f\n",
    median(boot_risk_mc$ci_width_pred, na.rm = TRUE),
    mean(boot_risk_mc$ci_width_pred,   na.rm = TRUE)))
} else {
  cat("Bootstrap could not run.\n")
}


# ─────────────────────────────────────────────────────────────────

if (!is.null(hold_probs) && exists("uncertainty_df")) {
  risk_df <- uncertainty_df %>%
    bind_cols(as.data.frame(avg_stack)) %>%
    mutate(
      risk_band = case_when(
        confidence >= 0.70 ~ "High confidence",
        confidence >= 0.50 ~ "Moderate confidence",
        TRUE               ~ "Low confidence — review recommended"
      )
    ) %>%
    select(sample_idx, true_subtype, pred_subtype, correct,
           confidence, disagreement, entropy, risk_band,
           all_of(lvls_hold))

  cat("\nRisk band distribution:\n"); print(table(risk_df$risk_band))

  knitr::kable(
    risk_df %>%
      arrange(desc(confidence)) %>%
      mutate(across(all_of(lvls_hold), ~round(.x, 3)),
             across(c(confidence, disagreement, entropy), round, 3)) %>%
      head(20),
    caption = "PAM50 subtype risk calculator — top 20 most confident predictions"
  )

  ggplot(risk_df, aes(x = pred_subtype, y = confidence, fill = pred_subtype)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 0.8) +
    theme_bw() +
    labs(x = "Predicted subtype", y = "Confidence (max P)",
         title = "Prediction confidence by PAM50 subtype") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))

  risk_df %>%
    arrange(desc(confidence)) %>%
    slice_head(n = 30) %>%
    mutate(sample_idx = factor(sample_idx, levels = sample_idx)) %>%
    pivot_longer(cols = all_of(lvls_hold),
                 names_to = "subtype", values_to = "prob") %>%
    ggplot(aes(x = sample_idx, y = prob, fill = subtype)) +
    geom_col(width = 0.9) +
    theme_bw(base_size = 9) +
    labs(x = "Sample (top 30 by confidence)", y = "P(subtype)", fill = "Subtype",
         title = "Predicted subtype probability profiles — top 30 most confident") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}


# ─────────────────────────────────────────────────────────────────

if (exists("risk_df")) {
  cat("=== CLASSIFIER SUMMARY ===\n\n")
  cat("Best model (LODO CV, ranked by mean MCC):", best_model_mc, "\n")
  cat(sprintf("  High-confidence:   %d / %d (%.0f%%)\n",
    sum(risk_df$risk_band == "High confidence"), nrow(risk_df),
    100 * mean(risk_df$risk_band == "High confidence", na.rm = TRUE)))
  cat(sprintf("  Moderate:          %d / %d (%.0f%%)\n",
    sum(risk_df$risk_band == "Moderate confidence"), nrow(risk_df),
    100 * mean(risk_df$risk_band == "Moderate confidence", na.rm = TRUE)))
  cat(sprintf("  Low (review):      %d / %d (%.0f%%)\n",
    sum(risk_df$risk_band == "Low confidence — review recommended"), nrow(risk_df),
    100 * mean(risk_df$risk_band == "Low confidence — review recommended",
               na.rm = TRUE)))
  cat(sprintf("\n  Median predicted confidence: %.1f%%\n",
    median(risk_df$confidence, na.rm = TRUE) * 100))
  if (exists("boot_arr_mc") && !is.null(boot_arr_mc)) {
    cat(sprintf("  Bootstrap 95%% CI (median width): %.1f%%\n",
      median(boot_risk_mc$ci_width_pred, na.rm = TRUE) * 100))
  }
}


# ─────────────────────────────────────────────────────────────────


# dir.create("models/multiclass", recursive=TRUE, showWarnings=FALSE)
# 
# ok_all    <- !is.na(p_merged_mc$label_mc) & p_merged_mc$label_mc!=""
# y_all     <- factor(p_merged_mc$label_mc[ok_all])
# e_all     <- e_merged_mc[,ok_all,drop=FALSE]
# batch_all <- p_merged_mc$dataset_id[ok_all]
# mod_all   <- model.matrix(~y_all)
# 
# grand_mean_all  <- rowMeans(e_all)
# batch_train_f   <- factor(as.character(batch_all))
# batch_means_all <- sapply(levels(batch_train_f), function(b) {
#   idx <- which(batch_train_f==b)
#   rowMeans(e_all[,idx,drop=FALSE])
# })
# 
# e_all_bc <- tryCatch(
#   sva::ComBat(dat=e_all,batch=batch_all,mod=mod_all,
#               par.prior=TRUE,mean.only=FALSE),
#   error=function(e) {
#     message("ComBat failed, using removeBatchEffect")
#     limma::removeBatchEffect(e_all,batch=batch_all,design=mod_all)
#   })
# 
# panel_genes <- select_genes_ovr(e_all_bc, y=y_all, n_per_class=50)
# cat("Final panel genes:", length(panel_genes), "\n")
# 
# X_all_bc  <- t(e_all_bc[panel_genes,,drop=FALSE])
# imp_final <- impute_fit(X_all_bc)
# X_all_i   <- impute_apply(X_all_bc,imp_final)
# sc_final  <- std_fit(X_all_i)
# X_all_s   <- std_apply(X_all_i,sc_final)
# 
# lvls_all <- levels(y_all); n_cl_all <- length(lvls_all)
# 
# # ── Fit and save the best model dynamically ───────────────────────────────────
# cat("Fitting and saving best model:", best_model_mc, "\n")
# 
# 
#  if (grepl("ElasticNet", best_model_mc)) {
#   fit_best <- glmnet::cv.glmnet(as.matrix(X_all_s), y_all,
#     family="multinomial", alpha=0.5, standardize=FALSE,
#     type.multinomial="grouped")
#   saveRDS(fit_best, "models/multiclass/fit_best_model.rds")
# 
# } else if (grepl("Lasso", best_model_mc)) {
#   fit_best <- glmnet::cv.glmnet(as.matrix(X_all_s), y_all,
#     family="multinomial", alpha=1, standardize=FALSE,
#     type.multinomial="grouped")
#   saveRDS(fit_best, "models/multiclass/fit_best_model.rds")
# 
# } else if (grepl("Ridge", best_model_mc)) {
#   fit_best <- glmnet::cv.glmnet(as.matrix(X_all_s), y_all,
#     family="multinomial", alpha=0, standardize=FALSE,
#     type.multinomial="grouped")
#   saveRDS(fit_best, "models/multiclass/fit_best_model.rds")
# 
# } else if (grepl("SVM_RBF", best_model_mc)) {
#   fit_best <- e1071::svm(X_all_s, y_all, kernel="radial", cost=1,
#     gamma=1/ncol(X_all_s), probability=TRUE)
#   saveRDS(fit_best, "models/multiclass/fit_best_model.rds")
# 
# } else if (grepl("SVM_Linear", best_model_mc)) {
#   fit_best <- e1071::svm(X_all_s, y_all, kernel="linear",
#     cost=1, probability=TRUE)
#   saveRDS(fit_best, "models/multiclass/fit_best_model.rds")
# 
# } else if (grepl("LDA", best_model_mc)) {
#   fit_best <- MASS::lda(x=X_all_s, grouping=y_all)
#   saveRDS(fit_best, "models/multiclass/fit_best_model.rds")
# 
# } else if (grepl("NaiveBayes", best_model_mc)) {
#   fit_best <- e1071::naiveBayes(x=X_all_s, y=y_all)
#   saveRDS(fit_best, "models/multiclass/fit_best_model.rds")
# 
# } else if (grepl("Multinom", best_model_mc)) {
#   fit_best <- nnet::multinom(y_all~., data=data.frame(y_all,X_all_s),
#     trace=FALSE, MaxNWts=100000)
#   saveRDS(fit_best, "models/multiclass/fit_best_model.rds")
# 
# } else {
#   warning("No save handler for model: ", best_model_mc,
#           ". Add a branch above.")
# }
# 
# # Save model identity so prediction function knows what was shipped
# saveRDS(best_model_mc, "models/multiclass/best_model_name.rds")
# 
# # Save preprocessing artefacts
# saveRDS(panel_genes, "models/multiclass/panel_genes.rds")
# saveRDS(list(grand_mean=grand_mean_all, batch_means=batch_means_all,
#              datasets_seen=levels(batch_train_f)),
#         "models/multiclass/batch_correction_params.rds")
# saveRDS(imp_final, "models/multiclass/imputation_params.rds")
# saveRDS(sc_final,  "models/multiclass/scale_params.rds")
# 
# cat("\nSaved to models/multiclass/:\n")
# cat(paste(" -", list.files("models/multiclass"), collapse="\n"), "\n")
# 
# if (exists("boot_risk_mc") && !is.null(boot_risk_mc))
#   saveRDS(boot_risk_mc, "models/multiclass/bootstrap_risk_scores.rds")
# if (exists("uncertainty_df") && !is.null(uncertainty_df))
#   saveRDS(uncertainty_df, "models/multiclass/model_disagreement.rds")
# 
# 
# # ─────────────────────────────────────────────────────────────────
# 
# # ── Lambda optimisation visualisation (glmnet models only) ───────────────────
# if (grepl("ElasticNet|Lasso|Ridge", best_model_mc)) {
#   fit_best <- readRDS("models/multiclass/fit_best_model.rds")
# 
#   lambda_df <- data.frame(
#     lambda     = fit_best$lambda,
#     cvm        = fit_best$cvm,
#     cvsd       = fit_best$cvsd,
#     log_lambda = log(fit_best$lambda)
#   ) %>%
#     mutate(lo = cvm - cvsd, hi = cvm + cvsd)
# 
#   lambda_min <- fit_best$lambda.min
#   lambda_1se <- fit_best$lambda.1se
# 
#   print(ggplot(lambda_df, aes(x = log_lambda, y = cvm)) +
#     geom_ribbon(aes(ymin = lo, ymax = hi),
#                 fill = "#2166AC", alpha = 0.15) +
#     geom_line(colour = "#2166AC", linewidth = 1.0) +
#     geom_vline(xintercept = log(lambda_min),
#                linetype = "dashed", colour = "#D7191C", linewidth = 0.8) +
#     annotate("text", x = log(lambda_min), y = max(lambda_df$hi),
#              label = sprintf("lambda.min\n(%.4f)", lambda_min),
#              colour = "#D7191C", hjust = -0.1, size = 3.5) +
#     theme_bw(base_size = 12) +
#     labs(
#       title    = paste0(best_model_mc, " — lambda optimisation via 10-fold CV"),
#       x = "log(lambda)", y = "Multinomial deviance (CV mean)"
#     ))
# 
#   nzero_df <- data.frame(
#     log_lambda = log(fit_best$lambda),
#     nzero      = fit_best$nzero
#   )
# 
#   print(ggplot(nzero_df, aes(x = log_lambda, y = nzero)) +
#     geom_line(colour = "#35964F", linewidth = 1.0) +
#     geom_point(size = 1.5, colour = "#35964F") +
#     geom_vline(xintercept = log(lambda_min),
#                linetype = "dashed", colour = "#D7191C", linewidth = 0.8) +
#     annotate("text", x = log(lambda_min),
#              y = max(nzero_df$nzero) * 0.9,
#              label = sprintf("lambda.min\n%d non-zero coefs",
#                nzero_df$nzero[which.min(abs(
#                  nzero_df$log_lambda - log(lambda_min)))]),
#              colour = "#D7191C", hjust = -0.1, size = 3.5) +
#     theme_bw(base_size = 12) +
#     labs(
#       title    = "Model sparsity along the lambda path",
#       subtitle = "Non-zero coefficients summed across all subtype contrasts",
#       x = "log(lambda)", y = "Number of non-zero coefficients"
#     ))
# } else {
#   cat("Lambda visualisation only available for glmnet models —",
#       best_model_mc, "was selected.\n")
# }


# ─────────────────────────────────────────────────────────────────

if (grepl("ElasticNet|Lasso|Ridge", best_model_mc)) {
  fit_final  <- readRDS("models/multiclass/fit_best_model.rds")
  coef_list  <- coef(fit_final, s = "lambda.min")

  coef_df <- bind_rows(lapply(names(coef_list), function(k) {
    cf <- as.matrix(coef_list[[k]])
    df <- data.frame(gene = rownames(cf), subtype = k,
                     coef = as.numeric(cf[, 1]), stringsAsFactors = FALSE)
    df <- df[df$gene != "(Intercept)", ]
    if (grepl("Ridge", best_model_mc))
      df %>% arrange(desc(abs(coef))) %>% slice_head(n = 15)
    else df[df$coef != 0, ]
  }))

  cat("Non-zero genes per subtype:\n"); print(table(coef_df$subtype))
  cat("Total unique genes:", length(unique(coef_df$gene)), "\n")

  coef_df %>%
    group_by(subtype) %>%
    slice_max(abs(coef), n = 10, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(direction = ifelse(coef > 0, "Positive", "Negative")) %>%
    ggplot(aes(x = reorder(gene, abs(coef)), y = coef, fill = direction)) +
    geom_col(alpha = 0.85) + coord_flip() +
    facet_wrap(~subtype, scales = "free_y") +
    scale_fill_manual(values = c("Positive" = "#2166AC",
                                 "Negative" = "#D7191C")) +
    theme_bw(base_size = 9) +
    labs(x = NULL, y = "Coefficient", fill = NULL,
         title = paste0("Top 10 genes per subtype — ", best_model_mc),
         subtitle = "Positive = upregulated in this subtype vs others") +
    theme(legend.position = "bottom")
}


# ─────────────────────────────────────────────────────────────────


# ── All metrics helper ─────────────────────────────────────────────────────────
compute_all_metrics <- function(truth, prob_mat) {
  lvls     <- colnames(prob_mat)
  truth_f  <- factor(as.character(truth), levels = lvls)
  pred_chr <- lvls[max.col(prob_mat)]
  pred_f   <- factor(pred_chr, levels = lvls)

  per_class <- lapply(lvls, function(k) {
    tp <- sum(pred_chr == k & truth == k, na.rm = TRUE)
    fp <- sum(pred_chr == k & truth != k, na.rm = TRUE)
    fn <- sum(pred_chr != k & truth == k, na.rm = TRUE)
    tn <- sum(pred_chr != k & truth != k, na.rm = TRUE)
    list(
      sens = if ((tp + fn) == 0) NA_real_ else tp / (tp + fn),
      spec = if ((tn + fp) == 0) NA_real_ else tn / (tn + fp),
      prec = if ((tp + fp) == 0) NA_real_ else tp / (tp + fp),
      npv  = if ((tn + fn) == 0) NA_real_ else tn / (tn + fn)
    )
  })

  y_bin <- sapply(lvls, function(k) as.integer(as.character(truth) == k))
  brier <- mean(rowSums((prob_mat - y_bin)^2))

  cm      <- table(Actual = truth_f, Predicted = pred_f)
  N       <- sum(cm)
  mcc_num <- N * sum(diag(cm)) - sum(rowSums(cm) * colSums(cm))
  mcc_den <- sqrt((N^2 - sum(colSums(cm)^2)) * (N^2 - sum(rowSums(cm)^2)))
  mcc     <- if (mcc_den == 0) 0 else mcc_num / mcc_den

  tibble::tibble(
    AUC            = round(macro_auc_ovr(truth, prob_mat), 3),
    Accuracy       = round(mean(pred_chr == truth, na.rm = TRUE), 3),
    `Balanced Acc` = round(balanced_accuracy(truth_f, pred_f), 3),
    Sensitivity    = round(mean(sapply(per_class, `[[`, "sens"), na.rm = TRUE), 3),
    Specificity    = round(mean(sapply(per_class, `[[`, "spec"), na.rm = TRUE), 3),
    Precision      = round(mean(sapply(per_class, `[[`, "prec"), na.rm = TRUE), 3),
    NPV            = round(mean(sapply(per_class, `[[`, "npv"),  na.rm = TRUE), 3),
    `F1 score`     = round(macro_f1(truth, pred_chr), 3),
    MCC            = round(mcc, 3),
    `Brier score`  = round(brier, 3)
  )
}

# ── Bootstrap hold-out metrics ─────────────────────────────────────────────────
bootstrap_holdout_metrics <- function(truth, prob_mat, B = 100) {
  lvls      <- colnames(prob_mat)
  strat_idx <- lapply(lvls, function(k) which(as.character(truth) == k))
  bind_rows(lapply(seq_len(B), function(b) {
    idx <- unlist(lapply(strat_idx, function(i) {
      if (length(i) == 0) return(integer(0))
      sample(i, length(i), replace = TRUE)
    }))
    m <- compute_all_metrics(truth[idx], prob_mat[idx, , drop = FALSE])
    m$iter <- b; m
  }))
}

# ── Forest plot ────────────────────────────────────────────────────────────────
forest_plot_metrics <- function(boot_df, point_metrics, n_per_class,
                                n_genes, B = 100) {
  metric_order <- c("Brier score", "MCC", "F1 score", "NPV", "Precision",
                    "Specificity", "Sensitivity", "Balanced Acc",
                    "Accuracy", "AUC")

  ci_df <- boot_df %>%
    pivot_longer(-iter, names_to = "metric", values_to = "value") %>%
    group_by(metric) %>%
    summarise(lo = quantile(value, 0.025, na.rm = TRUE),
              hi = quantile(value, 0.975, na.rm = TRUE), .groups = "drop")

  pt_df <- point_metrics %>%
    pivot_longer(everything(), names_to = "metric", values_to = "point")

  plot_df <- left_join(pt_df, ci_df, by = "metric") %>%
    filter(metric %in% metric_order) %>%
    mutate(metric = factor(metric, levels = metric_order))

  ggplot(plot_df, aes(x = point, y = metric)) +
    geom_vline(xintercept = c(0.25, 0.50, 0.75, 1.00),
               colour = "grey88", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi),
                   height = 0.25, colour = "#2166AC", linewidth = 1.0) +
    geom_point(size = 3.5, colour = "#D7191C") +
    scale_x_continuous(limits = c(0, 1.05),
                       breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          axis.title.y       = element_blank()) +
    labs(
      title    = paste0("Elastic Net Metric Uncertainty — ",
                        n_per_class, " genes/class (", n_genes, " total)"),
      subtitle = paste0(B, " class-stratified bootstrap resamples of the ",
                        "hold-out set (95% percentile interval)"),
      x = "Metric value"
    )
}

# ── Sweep ──────────────────────────────────────────────────────────────────────
n_grid         <- c(5, 10, 15, 20, 25, 30, 40, 50)
sweep_cache_v2 <- "../RDS_Files/multiclass/gene_sweep_correctedv4.rds"

set.seed(3888)
X_sw_dev      <- X_dev;  y_sw_dev      <- y_dev;  batch_sw_dev  <- batch_dev
X_sw_hold     <- X_hold; y_sw_hold     <- y_hold; batch_sw_hold <- batch_hold
cat("Dev:", nrow(X_sw_dev), "| Hold-out:", nrow(X_sw_hold), "\n")

if (file.exists(sweep_cache_v2)) {
  cat("Loading corrected sweep cache...\n")
  sweep_v2 <- readRDS(sweep_cache_v2)
} else {
  sweep_v2 <- lapply(n_grid, function(npc) {
    cat("\n===== n_per_class =", npc, "=====\n")

    # Panel gene count — single selection on full batch-corrected data, no leakage
    n_panel <- length(select_genes_ovr(e_mc_bc, y = y_mc, n_per_class = npc))
    cat("  Panel genes (full data):", n_panel, "\n")

    # LODO on full gene set — selection happens inside each fold
    res_npc <- run_lodo_multiclass(X_mc, y_mc, batch_mc,
                                   n_per_class = npc, use_weights = TRUE)
    if (nrow(res_npc$perf) == 0) { cat("  No LODO results\n"); return(NULL) }

    ba_npc <- res_npc$probs %>%
      group_by(fold, model) %>%
      group_modify(function(df, ...) {
        prob_cols <- setdiff(colnames(df), c("fold", "model", "truth"))
        lv <- prob_cols[prob_cols %in% levels(factor(df$truth))]
        if (length(lv) < 2) return(tibble::tibble(ba = NA_real_))
        pred <- lv[max.col(df[, lv, drop = FALSE])]
        tibble::tibble(ba = balanced_accuracy(df$truth, pred))
      }) %>% ungroup()

    ba_sum_npc <- ba_npc %>%
      group_by(model) %>%
      summarise(ba_mean = mean(ba, na.rm = TRUE), .groups = "drop")  # ← mean

    # Scorecard: mean across folds, sort by MCC then AUC
    sc_npc <- res_npc$perf %>%
      group_by(model) %>%
      summarise(
        m_auc = mean(macro_auc, na.rm = TRUE),   # ← mean
        m_f1  = mean(macro_f1,  na.rm = TRUE),   # ← mean
        m_mcc = mean(mcc,       na.rm = TRUE),   # ← mean MCC
        m_acc = mean(accuracy,  na.rm = TRUE),   # ← mean
        .groups = "drop"
      ) %>%
      left_join(ba_sum_npc, by = "model") %>%
      mutate(composite = (m_mcc + ba_mean + m_f1) / 3) %>%  # ← MCC in composite
      arrange(desc(m_mcc), desc(m_auc))  # ← sort by MCC

    target_model <- if ("ElasticNetMC"   %in% sc_npc$model) "ElasticNetMC"   else
                    if ("ElasticNetMC_w" %in% sc_npc$model) "ElasticNetMC_w" else
                    sc_npc$model[1]
    cat("  Reported model:", target_model, "\n")

    # LODO pooled metrics
    lodo_row     <- sc_npc %>% filter(model == target_model)
    pooled_probs <- res_npc$probs %>% filter(model == target_model)
    prob_cols    <- setdiff(colnames(pooled_probs), c("fold", "model", "truth"))
    lvls_pooled  <- prob_cols[prob_cols %in% unique(pooled_probs$truth)]
    pm_pooled    <- as.matrix(pooled_probs[, lvls_pooled, drop = FALSE])
    lodo_pooled  <- compute_all_metrics(pooled_probs$truth, pm_pooled)

    lodo_metrics <- lodo_pooled %>%
      mutate(AUC            = round(lodo_row$m_auc,  3),
             `F1 score`     = round(lodo_row$m_f1,   3),
             `Balanced Acc` = round(lodo_row$ba_mean, 3),
             MCC            = round(lodo_row$m_mcc,  3))

    cat("  LODO AUC:", lodo_metrics$AUC, "| MCC:", lodo_metrics$MCC, "\n")

    # Fit ElasticNet on dev, predict hold-out
    cb_hold <- combat_fit_apply(
      t(X_sw_dev), batch_sw_dev, model.matrix(~ y_sw_dev),
      t(X_sw_hold), batch_sw_hold
    )
    X_tr_bc <- t(cb_hold$train); X_te_bc <- t(cb_hold$test)

    gs_dev <- select_genes_ovr(t(X_tr_bc), y = y_sw_dev, n_per_class = npc)
    gs_dev <- intersect(gs_dev, colnames(X_tr_bc))
    if (length(gs_dev) < 10) { cat("  Too few genes\n"); return(NULL) }

    X_tr2 <- X_tr_bc[, gs_dev]; X_te2 <- X_te_bc[, gs_dev]
    imp2  <- impute_fit(X_tr2)
    X_tr2i <- impute_apply(X_tr2, imp2); X_te2i <- impute_apply(X_te2, imp2)
    sc2    <- std_fit(X_tr2i)
    X_tr2s <- std_apply(X_tr2i, sc2);   X_te2s <- std_apply(X_te2i, sc2)

    lvls_dev  <- levels(y_sw_dev); n_cl_dev <- length(lvls_dev)
    tab_dev   <- table(y_sw_dev);  n_tot_dev <- sum(tab_dev)
    sw_dev    <- setNames(vapply(lvls_dev, function(k)
      n_tot_dev / (n_cl_dev * max(1, tab_dev[k])), numeric(1)),
      lvls_dev)[as.character(y_sw_dev)]

    fit_hold <- tryCatch(
      glmnet::cv.glmnet(as.matrix(X_tr2s), y_sw_dev,
                        family = "multinomial", alpha = 0.5,
                        weights = sw_dev, standardize = FALSE,
                        type.multinomial = "grouped"),
      error = function(e) NULL
    )
    if (is.null(fit_hold)) { cat("  glmnet failed\n"); return(NULL) }

    pr_hold <- predict(fit_hold, as.matrix(X_te2s),
                       type = "response", s = "lambda.min")[,, 1]
    if (is.null(dim(pr_hold))) pr_hold <- t(as.matrix(pr_hold))
    lv_h    <- intersect(lvls_dev, colnames(pr_hold))
    pm_hold <- pr_hold[, lv_h, drop = FALSE]

    hold_metrics <- compute_all_metrics(as.character(y_sw_hold), pm_hold)
    cat("  Hold-out AUC:", hold_metrics$AUC, "| MCC:", hold_metrics$MCC, "\n")

    # Bootstrap hold-out metrics (B=100, test rows resampled, model fixed)
    cat("  Bootstrapping hold-out metrics (B=100)...\n")
    boot_metrics_df <- bootstrap_holdout_metrics(
      truth    = as.character(y_sw_hold),
      prob_mat = pm_hold, B = 100
    )

    # Confidence bands
    conf_hold  <- apply(pm_hold, 1, max)
    risk_bands <- dplyr::case_when(
      conf_hold >= 0.70 ~ "High confidence",
      conf_hold >= 0.50 ~ "Moderate confidence",
      TRUE              ~ "Low confidence"
    )

    # Model bootstrap for CI width (B=200)
    cat("  Model bootstrap for CI width (B=200)...\n")
    lambda_fix <- fit_hold$lambda.min
    n_tr <- nrow(X_tr2s); n_te <- nrow(X_te2s)

    ci_arr <- array(NA_real_, dim = c(n_te, n_cl_dev, 200),
                    dimnames = list(NULL, lvls_dev, NULL))
    for (b in seq_len(200)) {
      idx_b <- sample(n_tr, n_tr, replace = TRUE)
      y_b   <- droplevels(y_sw_dev[idx_b])
      X_b   <- X_tr2s[idx_b, , drop = FALSE]
      if (nlevels(y_b) < 2) next
      tab_b <- table(y_b); lvls_b <- levels(y_b)
      sw_b  <- setNames(vapply(lvls_b, function(k)
        sum(tab_b) / (length(lvls_b) * max(1, tab_b[k])), numeric(1)),
        lvls_b)[as.character(y_b)]
      lf <- tryCatch(
        glmnet::glmnet(as.matrix(X_b), y_b, family = "multinomial",
                       alpha = 0.5, weights = sw_b, standardize = FALSE,
                       lambda = lambda_fix, type.multinomial = "grouped"),
        error = function(e) NULL
      )
      if (is.null(lf)) next
      pr_b <- predict(lf, as.matrix(X_te2s), type = "response",
                      s = lambda_fix)[,, 1]
      if (is.null(dim(pr_b))) pr_b <- t(as.matrix(pr_b))
      lv_c <- intersect(lvls_dev, colnames(pr_b))
      ci_arr[, lv_c, b] <- pr_b[, lv_c]
    }

    pred_chr_hold <- lv_h[max.col(pm_hold)]
    ci_width_mat  <- apply(ci_arr, c(1, 2), function(x)
      diff(quantile(x, c(0.025, 0.975), na.rm = TRUE)))
    pred_idx <- match(pred_chr_hold, lvls_dev)
    ci_pred  <- vapply(seq_len(n_te), function(i)
      ci_width_mat[i, pred_idx[i]], numeric(1))

    unc_bands <- dplyr::case_when(
      ci_pred <= 0.10 ~ "Narrow",
      ci_pred <= 0.20 ~ "Moderate",
      TRUE            ~ "Wide"
    )
    cat("  Median CI width:", round(median(ci_pred, na.rm = TRUE), 3), "\n")

    list(
      n_per_class     = npc,
      n_panel_genes   = n_panel,
      best_model      = target_model,
      scorecard       = sc_npc,
      lodo_metrics    = lodo_metrics,
      hold_metrics    = hold_metrics,
      boot_metrics_df = boot_metrics_df,
      risk_bands      = table(risk_bands),
      unc_bands       = table(unc_bands),
      ci_median       = median(ci_pred, na.rm = TRUE)
    )
  })
  names(sweep_v2) <- paste0("n", n_grid)
  saveRDS(sweep_v2, sweep_cache_v2)
  cat("\nSaved to", sweep_cache_v2, "\n")
}

# ── Forest plots — one per n ───────────────────────────────────────────────────
for (r in sweep_v2) {
  if (is.null(r)) next
  print(forest_plot_metrics(
    boot_df       = r$boot_metrics_df,
    point_metrics = r$hold_metrics,
    n_per_class   = r$n_per_class,
    n_genes       = r$n_panel_genes,
    B             = 100
  ))
}

# ── Heatmap ────────────────────────────────────────────────────────────────────
metric_order <- c("Brier score", "MCC", "F1 score", "NPV", "Precision",
                  "Specificity", "Sensitivity", "Balanced Acc",
                  "Accuracy", "AUC")
metric_groups <- tibble::tibble(
  metric = metric_order,
  group  = c("Calibration", "Composite", "Composite", "Per-class",
             "Per-class",   "Per-class", "Per-class",
             "Discrimination", "Discrimination", "Discrimination")
)

heatmap_df <- bind_rows(lapply(sweep_v2, function(r) {
  if (is.null(r)) return(NULL)
  r$lodo_metrics %>%
    pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
    mutate(n_per_class   = r$n_per_class,
           n_panel_genes = r$n_panel_genes,
           panel_label   = paste0(r$n_per_class, "/class\n(",
                                  r$n_panel_genes, " genes)"))
})) %>%
  left_join(metric_groups, by = "metric") %>%
  filter(!is.na(group)) %>%
  mutate(metric        = factor(metric, levels = rev(metric_order)),
         value_display = ifelse(metric == "Brier score", 1 - value, value),
         label_val     = sprintf("%.3f", value))

ggplot(heatmap_df, aes(x = panel_label, y = metric, fill = value_display)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = label_val), size = 3, fontface = "bold") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradientn(
    colours = c("#D7191C", "#FDAE61", "#FFFFBF", "#A6D96A", "#1A9641"),
    limits  = c(0, 1), name = "Score\n(Brier inverted)") +
  scale_x_discrete(position = "top") +
  theme_bw(base_size = 11) +
  labs(title    = "All metrics across gene panel sizes — LODO pooled ElasticNet",
       subtitle = "Brier shown as 1 − Brier so green = better for every metric",
       x = NULL, y = NULL) +
  theme(strip.text.y     = element_text(angle = 0, face = "bold"),
        axis.text.x      = element_text(size = 9),
        panel.grid       = element_blank(),
        legend.position  = "right")

# ── Line traces ────────────────────────────────────────────────────────────────
line_df <- bind_rows(lapply(sweep_v2, function(r) {
  if (is.null(r)) return(NULL)
  bind_rows(
    r$lodo_metrics %>%
      pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
      mutate(source = "LODO"),
    r$hold_metrics %>%
      pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
      mutate(source = "Hold-out")
  ) %>% mutate(n_panel_genes = r$n_panel_genes, n_per_class = r$n_per_class)
})) %>%
  left_join(metric_groups, by = "metric") %>%
  filter(!is.na(group), metric != "Brier score")

ggplot(line_df, aes(x = n_panel_genes, y = value,
                    colour = source, linetype = source,
                    group = interaction(metric, source))) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  facet_wrap(~ metric, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = c("LODO" = "#2166AC", "Hold-out" = "#E87722")) +
  scale_linetype_manual(values = c("LODO" = "solid", "Hold-out" = "dashed")) +
  theme_bw(base_size = 10) +
  labs(title    = "LODO vs Hold-out across panel sizes",
       subtitle = "Large gap = optimism in LODO estimate",
       x = "Total unique panel genes", y = "Score",
       colour = NULL, linetype = NULL) +
  theme(legend.position = "bottom")

# ── Summary table ──────────────────────────────────────────────────────────────
sweep_summary_v2 <- bind_rows(lapply(sweep_v2, function(r) {
  if (is.null(r)) return(NULL)
  tibble::tibble(
    `N/class`     = r$n_per_class,
    `Panel genes` = r$n_panel_genes,
    `LODO AUC`    = r$lodo_metrics$AUC,
    `LODO MCC`    = r$lodo_metrics$MCC,
    `LODO BA`     = r$lodo_metrics$`Balanced Acc`,
    `LODO F1`     = r$lodo_metrics$`F1 score`,
    `LODO Brier`  = r$lodo_metrics$`Brier score`,
    `Hold AUC`    = r$hold_metrics$AUC,
    `Hold MCC`    = r$hold_metrics$MCC,
    `Hold BA`     = r$hold_metrics$`Balanced Acc`,
    `% High conf` = if (!is.null(r$risk_bands))
      round(100 * r$risk_bands["High confidence"] / sum(r$risk_bands), 1) else NA,
    `% Narrow CI` = if (!is.null(r$unc_bands))
      round(100 * r$unc_bands["Narrow"] / sum(r$unc_bands), 1) else NA,
    `Median CI`   = round(r$ci_median, 3)
  )
}))

knitr::kable(sweep_summary_v2,
  caption = "Corrected gene panel sweep — ElasticNet, no leakage, sorted by MCC")

# ── Plot 1: AUC flat, but MCC and BA still climb ──────────────────────────────
# This is the core argument: AUC saturates early, other metrics do not.
sweep_summary_v2 %>%
  select(`N/class`, `Panel genes`, `LODO AUC`, `LODO MCC`, `LODO BA`, `LODO F1`) %>%
  pivot_longer(cols = c(`LODO AUC`, `LODO MCC`, `LODO BA`, `LODO F1`),
               names_to = "metric", values_to = "value") %>%
  mutate(
    metric     = factor(metric,
                        levels = c("LODO AUC", "LODO MCC", "LODO BA", "LODO F1")),
    saturated  = metric == "LODO AUC"
  ) %>%
  ggplot(aes(x = `Panel genes`, y = value,
             colour = metric, group = metric)) +
  geom_line(aes(linewidth = saturated, linetype = saturated)) +
  geom_point(size = 3) +
  geom_text(data = ~filter(.x, `Panel genes` == max(`Panel genes`)),
            aes(label = sprintf("%.3f", value)),
            hjust = -0.2, size = 3.2) +
  scale_linewidth_manual(values = c("TRUE" = 0.6, "FALSE" = 1.2),
                         guide = "none") +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                        guide = "none") +
  scale_colour_manual(values = c(
    "LODO AUC" = "grey60",
    "LODO MCC" = "#D7191C",
    "LODO BA"  = "#2166AC",
    "LODO F1"  = "#35964F"
  )) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_y_continuous(limits = c(0, 1)) +
  annotate("text", x = min(sweep_summary_v2$`Panel genes`) + 5,
           y = 0.97, label = "AUC saturates early →",
           colour = "grey50", size = 3.2, hjust = 0, fontface = "italic") +
  theme_bw(base_size = 12) +
  labs(
    title    = "AUC saturates; MCC and Balanced Accuracy continue to improve",
    subtitle = "Dashed grey = AUC (ranking metric). Solid = threshold-sensitive metrics.",
    x = "Total unique panel genes", y = "Score", colour = NULL
  ) +
  theme(legend.position = "bottom")


# ── Plot 2: Calibration and uncertainty improve with panel size ───────────────
# Brier score and CI width are the calibration arguments for n=50.
sweep_summary_v2 %>%
  select(`N/class`, `Panel genes`, `LODO Brier`, `Median CI`) %>%
  pivot_longer(cols = c(`LODO Brier`, `Median CI`),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
    `LODO Brier` = "Brier score (lower = better)",
    `Median CI`  = "Median bootstrap 95% CI width (lower = better)"
  )) %>%
  ggplot(aes(x = `Panel genes`, y = value,
             colour = metric, group = metric)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_text(data = ~filter(.x, `Panel genes` == max(`Panel genes`)),
            aes(label = sprintf("%.3f", value)),
            hjust = -0.2, size = 3.2) +
  scale_colour_manual(values = c(
    "Brier score (lower = better)"                       = "#E87722",
    "Median bootstrap 95% CI width (lower = better)"    = "#6A3D9A"
  )) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_bw(base_size = 12) +
  labs(
    title    = "Calibration and prediction certainty improve with panel size",
    subtitle = "Both metrics are minimised — larger panel → better calibrated probabilities and narrower CIs",
    x = "Total unique panel genes", y = "Score", colour = NULL
  ) +
  theme(legend.position = "bottom")


# ── Plot 3: High-confidence % and narrow CI % ─────────────────────────────────
# Clinical utility argument: more genes = more patients get a confident call.
sweep_summary_v2 %>%
  select(`N/class`, `Panel genes`, `% High conf`, `% Narrow CI`) %>%
  pivot_longer(cols = c(`% High conf`, `% Narrow CI`),
               names_to = "metric", values_to = "pct") %>%
  mutate(metric = recode(metric,
    `% High conf`  = "% patients: high confidence (≥0.70)",
    `% Narrow CI`  = "% patients: narrow CI (width ≤0.10)"
  )) %>%
  ggplot(aes(x = `Panel genes`, y = pct,
             colour = metric, group = metric)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  geom_text(data = ~filter(.x, `Panel genes` == max(`Panel genes`)),
            aes(label = paste0(round(pct, 1), "%")),
            hjust = -0.2, size = 3.2) +
  scale_colour_manual(values = c(
    "% patients: high confidence (≥0.70)"  = "#2166AC",
    "% patients: narrow CI (width ≤0.10)"  = "#35964F"
  )) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_bw(base_size = 12) +
  labs(
    title    = "Clinical utility: more genes → more patients receive a confident call",
    subtitle = "Larger panel reduces the proportion flagged for manual review",
    x = "Total unique panel genes", y = "% of hold-out patients",
    colour = NULL
  ) +
  theme(legend.position = "bottom")


# ── Plot 4: Elbow detection — marginal gain per additional gene ───────────────
# Shows the gain is real early on but the n=50 plateau is defensible.
sweep_summary_v2 %>%
  arrange(`Panel genes`) %>%
  mutate(
    mcc_gain = c(NA, diff(`LODO MCC`)),
    ba_gain  = c(NA, diff(`LODO BA`)),
    genes_added = c(NA, diff(`Panel genes`)),
    mcc_per_gene = mcc_gain / genes_added,
    ba_per_gene  = ba_gain  / genes_added
  ) %>%
  filter(!is.na(mcc_gain)) %>%
  select(`Panel genes`, mcc_per_gene, ba_per_gene) %>%
  pivot_longer(cols = c(mcc_per_gene, ba_per_gene),
               names_to = "metric", values_to = "gain_per_gene") %>%
  mutate(metric = recode(metric,
    mcc_per_gene = "MCC gain per added gene",
    ba_per_gene  = "Balanced Acc gain per added gene"
  )) %>%
  ggplot(aes(x = `Panel genes`, y = gain_per_gene,
             colour = metric, group = metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 1.0) +
  geom_point(size = 3) +
  scale_colour_manual(values = c(
    "MCC gain per added gene"            = "#D7191C",
    "Balanced Acc gain per added gene"   = "#2166AC"
  )) +
  theme_bw(base_size = 12) +
  labs(
    title    = "Marginal gain per additional gene — diminishing but positive at n=50",
    subtitle = "Gains approach zero but remain positive, supporting n=50 over early stopping",
    x = "Total unique panel genes", y = "Metric gain per gene added",
    colour = NULL
  ) +
  theme(legend.position = "bottom")


# ═══════════════════════════════════════════════════════════════════════════════
# ADDITIONAL RDS SAVES FOR REPORT REPRODUCIBILITY
# These saves must run AFTER the scorecard chunk and AFTER the gene sweep chunk.
# ═══════════════════════════════════════════════════════════════════════════════

dir.create("../RDS_Files/multiclass", showWarnings = FALSE, recursive = TRUE)

# ── LODO performance and probs ────────────────────────────────────────────────
saveRDS(res_lodo$perf,  "../RDS_Files/multiclasslodo_perf_multiclass.rds")
saveRDS(res_lodo$probs, "../RDS_Files/multiclass/lodo_probs_multiclass.rds")

# ── Scorecard and BA summary ──────────────────────────────────────────────────
saveRDS(scorecard,   "../RDS_Files/multiclass/scorecard_multiclass.rds")
saveRDS(ba_summary,  "../RDS_Files/multiclass/ba_summary_multiclass.rds")

# ── Best model name (also in models/ but duplicate here for report path) ──────
saveRDS(best_model_mc, "../RDS_Files/multiclass/best_model_mc_multiclass.rds")

# ── Pooled predictions and per-class metrics ──────────────────────────────────
saveRDS(preds_best,   "../RDS_Files/multiclass/preds_best_multiclass.rds")
saveRDS(sens_spec_tbl,"../RDS_Files/multiclass/sens_spec_tbl_multiclass.rds")

# ── y_mc levels and hold-out labels ───────────────────────────────────────────
saveRDS(y_mc,     "../RDS_Files/multiclass/y_mc_multiclass.rds")
saveRDS(all_lvls, "../RDS_Files/multiclass/all_lvls_multiclass.rds")
saveRDS(y_hold,   "../RDS_Files/multiclass/y_hold_multiclass.rds")

# boot_risk_mc already saved to models/multiclass/ — copy to results too
if (exists("boot_risk_mc") && !is.null(boot_risk_mc))
  saveRDS(boot_risk_mc, "../RDS_Files/multiclass/boot_risk_mc_multiclass.rds")

cat("\nAll multiclass result RDS files saved.\n")
cat("Files in ../RDS_Files/multiclass/:\n")
cat(paste(" -", list.files("../RDS_Files/multiclass", pattern="\\.rds$")), sep="\n")
