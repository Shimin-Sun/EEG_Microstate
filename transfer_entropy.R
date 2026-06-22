if (!require(dplyr)) install.packages("dplyr")
if (!require(RTransferEntropy)) install.packages("RTransferEntropy")
if (!require(openxlsx)) install.packages("openxlsx")
if (!require(parallel)) install.packages("parallel")
if (!require(doParallel)) install.packages("doParallel")
if (!require(foreach)) install.packages("foreach")

library(dplyr)
library(RTransferEntropy)
library(openxlsx)
library(parallel)
library(doParallel)
library(foreach)

base_path <- "path/to/EC_sequences"
output_dir <- "path/to/output/transfer_entropy"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

hc_range <- 1:162
mdd_range <- 6:416
condition <- "EC"
microstates <- c("A", "B", "C", "D")
sampling_rate <- 250

band_combinations <- list(
  c("Broadband", "Delta"),
  c("Broadband", "Theta"),
  c("Broadband", "Alpha"),
  c("Broadband", "Beta"),
  c("Broadband", "Gamma")
)

calculate_cohens_d <- function(group1, group2) {
  n1 <- length(group1)
  n2 <- length(group2)
  sd1 <- sd(group1, na.rm = TRUE)
  sd2 <- sd(group2, na.rm = TRUE)
  pooled_var <- ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)
  pooled_sd <- sqrt(pooled_var)
  mean_diff <- mean(group2, na.rm = TRUE) - mean(group1, na.rm = TRUE)
  if (pooled_sd == 0) return(NA)
  return(mean_diff / pooled_sd)
}

validate_parameters <- function() {
  if (!dir.exists(base_path)) stop("Base path does not exist: ", base_path)
  if (length(microstates) == 0) stop("Microstate list is empty")
  if (length(hc_range) == 0 || length(mdd_range) == 0) stop("Subject range error")
  if (!condition %in% c("EC", "EO")) stop("Condition must be 'EC' or 'EO'")
  if (length(band_combinations) == 0) stop("Band combination list is empty")
  
  if (!dir.exists(output_dir)) {
    tryCatch({
      dir.create(output_dir, recursive = TRUE)
    }, error = function(e) stop("Cannot create output directory: ", e$message))
  }
  
  test_file <- file.path(output_dir, "test_write.txt")
  tryCatch({
    writeLines("test", test_file)
    file.remove(test_file)
  }, error = function(e) stop("No write permission: ", output_dir))
}

check_file_existence <- function(subject_id, group, base_path, condition, band1, band2) {
  tryCatch({
    band1_file <- file.path(base_path, band1, sprintf("%s%s_%s_%s.set_dynamics.txt", group, subject_id, condition, band1))
    band2_file <- file.path(base_path, band2, sprintf("%s%s_%s_%s.set_dynamics.txt", group, subject_id, condition, band2))
    if (!file.exists(band1_file) || !file.exists(band2_file)) return(FALSE)
    if (file.info(band1_file)$size == 0 || file.info(band2_file)$size == 0) return(FALSE)
    return(TRUE)
  }, error = function(e) return(FALSE))
}

read_microstate_file <- function(file_path) {
  file_path <- normalizePath(file_path, mustWork = FALSE)
  if (!file.exists(file_path)) return(NULL)
  tryCatch({
    data <- read.table(file_path, header = TRUE, sep = "", na.strings = c("NA", ""), stringsAsFactors = FALSE)
    if (nrow(data) == 0) return(NULL)
    required_cols <- c("Time", "A", "B", "C", "D")
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) > 0) return(NULL)
    return(data)
  }, error = function(e) return(NULL))
}

convert_to_labels <- function(data) {
  intensity_cols <- c("A", "B", "C", "D")
  labels <- apply(data[, intensity_cols], 1, function(row) {
    non_zero <- which(row != 0)
    if (length(non_zero) == 0) return(NA)
    else if (length(non_zero) == 1) return(intensity_cols[non_zero])
    else {
      abs_values <- abs(row)
      return(intensity_cols[which.max(abs_values)])
    }
  })
  return(data.frame(Time = data$Time, Label = labels, stringsAsFactors = FALSE))
}

convert_to_numeric_sequence <- function(labels_df) {
  state_mapping <- c("A" = 1, "B" = 2, "C" = 3, "D" = 4)
  numeric_seq <- sapply(labels_df$Label, function(label) {
    if (is.na(label)) return(0)
    else return(state_mapping[label])
  })
  return(numeric_seq)
}

perform_transfer_entropy_analysis <- function(seq1, seq2, lx = 1, ly = 1, nboot = 1000) {
  min_length <- min(length(seq1), length(seq2))
  if (min_length < 50) return(NULL)
  seq1_trimmed <- seq1[1:min_length]
  seq2_trimmed <- seq2[1:min_length]
  valid_indices <- !is.na(seq1_trimmed) & !is.na(seq2_trimmed) & seq1_trimmed != 0 & seq2_trimmed != 0
  seq1_valid <- seq1_trimmed[valid_indices]
  seq2_valid <- seq2_trimmed[valid_indices]
  if (length(seq1_valid) < 50) return(NULL)
  if (length(unique(seq1_valid)) < 2 || length(unique(seq2_valid)) < 2) return(NULL)
  
  tryCatch({
    te_result <- transfer_entropy(
      x = seq1_valid,
      y = seq2_valid,
      lx = lx,
      ly = ly,
      nboot = nboot,
      type = "limits",
      limits = c(0.5, 1.5, 2.5, 3.5, 4.5)
    )
    return(list(
      te_1to2 = te_result$coef[1, 1], te_2to1 = te_result$coef[2, 1],
      p_value_1to2 = te_result$coef[1, 4], p_value_2to1 = te_result$coef[2, 4],
      eff_te_1to2 = te_result$coef[1, 2], eff_te_2to1 = te_result$coef[2, 2],
      n_effective = length(seq1_valid)
    ))
  }, error = function(e) return(NULL))
}

process_single_subject_te <- function(subject_id, group, base_path, condition, band1, band2, microstates, lx = 1, ly = 1, nboot = 1000) {
  tryCatch({
    if (!check_file_existence(subject_id, group, base_path, condition, band1, band2)) return(NULL)
    
    band1_file <- file.path(base_path, band1, sprintf("%s%s_%s_%s.set_dynamics.txt", group, subject_id, condition, band1))
    band2_file <- file.path(base_path, band2, sprintf("%s%s_%s_%s.set_dynamics.txt", group, subject_id, condition, band2))
    
    data_band1 <- read_microstate_file(band1_file)
    data_band2 <- read_microstate_file(band2_file)
    if (is.null(data_band1) || is.null(data_band2)) return(NULL)
    if (!all.equal(data_band1$Time, data_band2$Time, tolerance = 1e-6)) return(NULL)
    
    labels_band1 <- convert_to_labels(data_band1)
    labels_band2 <- convert_to_labels(data_band2)
    seq_band1 <- convert_to_numeric_sequence(labels_band1)
    seq_band2 <- convert_to_numeric_sequence(labels_band2)
    
    te_result <- perform_transfer_entropy_analysis(seq_band1, seq_band2, lx, ly, nboot)
    if (is.null(te_result)) return(NULL)
    
    return(list(
      subject_id = paste0(group, subject_id), group = group, band1 = band1, band2 = band2,
      n_total_band1 = length(seq_band1), n_total_band2 = length(seq_band2), te_results = te_result
    ))
  }, error = function(e) return(NULL))
}

perform_band_specific_fdr <- function(combo_results) {
  all_p_values <- c()
  all_tests <- c()
  all_subjects <- c()
  
  for (subject_result in combo_results) {
    if (is.null(subject_result$te_results)) next
    p1 <- subject_result$te_results$p_value_1to2
    if (!is.na(p1)) {
      all_p_values <- c(all_p_values, p1)
      all_tests <- c(all_tests, paste0(subject_result$band1, "->", subject_result$band2))
      all_subjects <- c(all_subjects, subject_result$subject_id)
    }
    p2 <- subject_result$te_results$p_value_2to1
    if (!is.na(p2)) {
      all_p_values <- c(all_p_values, p2)
      all_tests <- c(all_tests, paste0(subject_result$band2, "->", subject_result$band1))
      all_subjects <- c(all_subjects, subject_result$subject_id)
    }
  }
  
  if (length(all_p_values) == 0) return(NULL)
  fdr_adjusted <- p.adjust(all_p_values, method = "fdr")
  return(data.frame(
    subject_id = all_subjects, test = all_tests, raw_p = all_p_values, fdr_adjusted = fdr_adjusted,
    significant = fdr_adjusted < 0.05, stringsAsFactors = FALSE
  ))
}

perform_group_comparison <- function(combo_results, combo_name) {
  if (length(combo_results) == 0) return(NULL)
  
  test_directions <- c(
    paste0(combo_results[[1]]$band1, "->", combo_results[[1]]$band2),
    paste0(combo_results[[1]]$band2, "->", combo_results[[1]]$band1)
  )
  
  group_comparison_results <- data.frame()
  
  for (test_dir in test_directions) {
    hc_te <- c(); hc_eff_te <- c(); hc_p <- c()
    mdd_te <- c(); mdd_eff_te <- c(); mdd_p <- c()
    
    for (subject_result in combo_results) {
      if (is.null(subject_result$te_results)) next
      if (test_dir == paste0(subject_result$band1, "->", subject_result$band2)) {
        te_val <- subject_result$te_results$te_1to2
        eff_te_val <- subject_result$te_results$eff_te_1to2
        p_val <- subject_result$te_results$p_value_1to2
      } else if (test_dir == paste0(subject_result$band2, "->", subject_result$band1)) {
        te_val <- subject_result$te_results$te_2to1
        eff_te_val <- subject_result$te_results$eff_te_2to1
        p_val <- subject_result$te_results$p_value_2to1
      } else { next }
      
      if (subject_result$group == "HC") {
        hc_te <- c(hc_te, te_val); hc_eff_te <- c(hc_eff_te, eff_te_val); hc_p <- c(hc_p, p_val)
      } else if (subject_result$group == "MDD") {
        mdd_te <- c(mdd_te, te_val); mdd_eff_te <- c(mdd_eff_te, eff_te_val); mdd_p <- c(mdd_p, p_val)
      }
    }
    
    if (length(hc_te) >= 3 && length(mdd_te) >= 3) {
      hc_te_mean <- mean(hc_te, na.rm = TRUE); hc_te_sd <- sd(hc_te, na.rm = TRUE)
      mdd_te_mean <- mean(mdd_te, na.rm = TRUE); mdd_te_sd <- sd(mdd_te, na.rm = TRUE)
      te_mean_diff <- mdd_te_mean - hc_te_mean
      var_test_te <- tryCatch(var.test(hc_te, mdd_te), error=function(e) list(p.value=NA))
      var_equal_te <- ifelse(is.na(var_test_te$p.value), FALSE, var_test_te$p.value > 0.05)
      te_t_test <- tryCatch(t.test(hc_te, mdd_te, var.equal=var_equal_te), error=function(e) list(statistic=NA,p.value=NA,stderr=NA))
      te_t <- as.numeric(-te_t_test$statistic); te_se <- te_t_test$stderr; te_p_raw <- te_t_test$p.value
      te_test_type <- ifelse(var_equal_te, "Student's t-test", "Welch's t-test")
      te_cohens_d <- calculate_cohens_d(hc_te, mdd_te)
      
      hc_eff_te_mean <- mean(hc_eff_te, na.rm = TRUE); hc_eff_te_sd <- sd(hc_eff_te, na.rm = TRUE)
      mdd_eff_te_mean <- mean(mdd_eff_te, na.rm = TRUE); mdd_eff_te_sd <- sd(mdd_eff_te, na.rm = TRUE)
      eff_te_mean_diff <- mdd_eff_te_mean - hc_eff_te_mean
      var_test_eff <- tryCatch(var.test(hc_eff_te, mdd_eff_te), error=function(e) list(p.value=NA))
      var_equal_eff <- ifelse(is.na(var_test_eff$p.value), FALSE, var_test_eff$p.value > 0.05)
      eff_te_t_test <- tryCatch(t.test(hc_eff_te, mdd_eff_te, var.equal=var_equal_eff), error=function(e) list(statistic=NA,p.value=NA,stderr=NA))
      eff_te_t <- as.numeric(-eff_te_t_test$statistic); eff_te_se <- eff_te_t_test$stderr; eff_te_p_raw <- eff_te_t_test$p.value
      eff_te_test_type <- ifelse(var_equal_eff, "Student's t-test", "Welch's t-test")
      eff_te_cohens_d <- calculate_cohens_d(hc_eff_te, mdd_eff_te)
      
      hc_sig_count <- sum(hc_p < 0.05, na.rm = TRUE); mdd_sig_count <- sum(mdd_p < 0.05, na.rm = TRUE)
      hc_total <- length(hc_p); mdd_total <- length(mdd_p)
      prop_test_result <- tryCatch({
        if (hc_sig_count > 0 || mdd_sig_count > 0) prop.test(x=c(hc_sig_count,mdd_sig_count),n=c(hc_total,mdd_total))
        else list(p.value=NA)
      }, error=function(e) list(p.value=NA))
      prop_p_raw <- prop_test_result$p.value
      
      group_comparison_results <- rbind(group_comparison_results,
                                        data.frame(
                                          test_direction = test_dir, band_combo = combo_name,
                                          hc_sample_size = length(hc_te), mdd_sample_size = length(mdd_te),
                                          hc_te_mean = hc_te_mean, hc_te_sd = hc_te_sd,
                                          mdd_te_mean = mdd_te_mean, mdd_te_sd = mdd_te_sd,
                                          te_mean_diff = te_mean_diff, te_se = te_se,
                                          te_t_value = te_t, te_test_type = te_test_type,
                                          te_p_raw = te_p_raw, te_p_fdr = NA,
                                          te_cohens_d = te_cohens_d, te_significant = FALSE,
                                          hc_eff_te_mean = hc_eff_te_mean, hc_eff_te_sd = hc_eff_te_sd,
                                          mdd_eff_te_mean = mdd_eff_te_mean, mdd_eff_te_sd = mdd_eff_te_sd,
                                          eff_te_mean_diff = eff_te_mean_diff, eff_te_se = eff_te_se,
                                          eff_te_t_value = eff_te_t, eff_te_test_type = eff_te_test_type,
                                          eff_te_p_raw = eff_te_p_raw, eff_te_p_fdr = NA,
                                          eff_te_cohens_d = eff_te_cohens_d, eff_te_significant = FALSE,
                                          hc_sig_prop = hc_sig_count/hc_total, mdd_sig_prop = mdd_sig_count/mdd_total,
                                          prop_p_raw = prop_p_raw, prop_p_fdr = NA, prop_significant = FALSE,
                                          stringsAsFactors = FALSE
                                        ))
    }
  }
  
  if (nrow(group_comparison_results) > 0) {
    group_comparison_results$te_p_fdr <- p.adjust(group_comparison_results$te_p_raw, "fdr")
    group_comparison_results$te_significant <- group_comparison_results$te_p_fdr < 0.05
    group_comparison_results$eff_te_p_fdr <- p.adjust(group_comparison_results$eff_te_p_raw, "fdr")
    group_comparison_results$eff_te_significant <- group_comparison_results$eff_te_p_fdr < 0.05
    group_comparison_results$prop_p_fdr <- p.adjust(group_comparison_results$prop_p_raw, "fdr")
    group_comparison_results$prop_significant <- group_comparison_results$prop_p_fdr < 0.05
  }
  return(group_comparison_results)
}

process_single_group <- function(band1, band2, group, subject_range, base_path, condition, microstates, lx, ly, nboot, combo_output_dir, combo_name) {
  log_file <- file.path(combo_output_dir, paste0(combo_name, "_", group, "_log.txt"))
  cat(paste0("=== ", combo_name, " ", group, " log ===\n", "Start: ", Sys.time(), "\n\n"), file = log_file)
  
  group_results <- list()
  
  for (i in seq_along(subject_range)) {
    subject_id <- subject_range[i]
    if (i %% 10 == 0) gc()
    
    result <- process_single_subject_te(
      subject_id = subject_id, group = group, base_path = base_path,
      condition = condition, band1 = band1, band2 = band2, microstates = microstates,
      lx = lx, ly = ly, nboot = nboot
    )
    
    if (!is.null(result)) {
      log_msg <- sprintf("[%s] SUCCESS: %s%s (n_eff=%d)\n", format(Sys.time(), "%H:%M:%S"), group, subject_id, result$te_results$n_effective)
      group_results[[length(group_results) + 1]] <- result
    } else {
      log_msg <- sprintf("[%s] FAIL: %s%s\n", format(Sys.time(), "%H:%M:%S"), group, subject_id)
    }
    cat(log_msg, file = log_file, append = TRUE)
    if (i %% 5 == 0) cat(sprintf("[%s-%s] Progress: %d/%d\n", group, combo_name, i, length(subject_range)))
  }
  
  cat(sprintf("[DONE] %s-%s valid subjects: %d\n", group, combo_name, length(group_results)))
  return(group_results)
}

finalize_band_combo <- function(band1, band2, hc_results, mdd_results, combo_output_dir, combo_name) {
  combo_results <- c(hc_results, mdd_results)
  if (length(combo_results) == 0) {
    message("[WARN] ", combo_name, " no valid results")
    return(NULL)
  }
  
  fdr_results <- perform_band_specific_fdr(combo_results)
  group_comp <- perform_group_comparison(combo_results, combo_name)
  
  if (!is.null(fdr_results)) write.csv(fdr_results, file.path(combo_output_dir, paste0(combo_name, "_FDR.csv")), row.names = F)
  if (!is.null(group_comp)) write.csv(group_comp, file.path(combo_output_dir, paste0(combo_name, "_group_diff.csv")), row.names = F)
  
  detailed_rows <- list()
  for (subj_res in combo_results) {
    detailed_rows[[length(detailed_rows) + 1]] <- data.frame(
      subject_id = subj_res$subject_id, group = subj_res$group,
      direction = paste0(subj_res$band1, "->", subj_res$band2),
      transfer_entropy = subj_res$te_results$te_1to2,
      effective_te = subj_res$te_results$eff_te_1to2,
      p_value = subj_res$te_results$p_value_1to2,
      n_effective = subj_res$te_results$n_effective, stringsAsFactors = F)
    detailed_rows[[length(detailed_rows) + 1]] <- data.frame(
      subject_id = subj_res$subject_id, group = subj_res$group,
      direction = paste0(subj_res$band2, "->", subj_res$band1),
      transfer_entropy = subj_res$te_results$te_2to1,
      effective_te = subj_res$te_results$eff_te_2to1,
      p_value = subj_res$te_results$p_value_2to1,
      n_effective = subj_res$te_results$n_effective, stringsAsFactors = F)
  }
  
  detailed_df <- do.call(rbind, detailed_rows)
  write.csv(detailed_df, file.path(combo_output_dir, paste0(combo_name, "_detailed.csv")), row.names = F)
  
  summary_df <- detailed_df %>% group_by(group, direction) %>%
    summarise(n = n(), te_mean = mean(transfer_entropy, na.rm = T), te_sd = sd(transfer_entropy, na.rm = T),
              eff_te_mean = mean(effective_te, na.rm = T), p_mean = mean(p_value, na.rm = T),
              sig_prop = mean(p_value < 0.05, na.rm = T), .groups = "drop")
  write.csv(summary_df, file.path(combo_output_dir, paste0(combo_name, "_group_summary.csv")), row.names = F)
  
  saveRDS(list(detailed = detailed_df, summary = summary_df, fdr = fdr_results, group_comp = group_comp),
          file.path(combo_output_dir, paste0(combo_name, "_results.rds")))
  
  return(data.frame(band_combo = combo_name, n_subjects = length(combo_results), stringsAsFactors = F))
}

validate_parameters()
start_time <- Sys.time()
lx <- 1; ly <- 1; nboot <- 1000

tasks <- list()
for (combo in band_combinations) {
  b1 <- combo[1]
  b2 <- combo[2]
  cname <- paste0(b1, "-", b2)
  
  combo_dir <- file.path(output_dir, cname)
  if (!dir.exists(combo_dir)) dir.create(combo_dir, recursive = TRUE)
  
  tasks[[length(tasks) + 1]] <- list(band1 = b1, band2 = b2, group = "HC", range = hc_range, name = cname, dir = combo_dir)
  tasks[[length(tasks) + 1]] <- list(band1 = b1, band2 = b2, group = "MDD", range = mdd_range, name = cname, dir = combo_dir)
}

cat("Starting parallel processing...\n")
cat("Total band combinations:", length(band_combinations), "\n")
cat("Total tasks:", length(tasks), "\n")

max_cores <- detectCores() - 1
n_cores <- min(max_cores, length(tasks))
cat("Using cores:", n_cores, "\n\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterExport(cl, c("base_path", "condition", "microstates", "lx", "ly", "nboot", "output_dir",
                    "check_file_existence", "read_microstate_file", "convert_to_labels",
                    "convert_to_numeric_sequence", "perform_transfer_entropy_analysis",
                    "process_single_subject_te", "process_single_group"),
              envir = environment())
clusterEvalQ(cl, { library(dplyr); library(RTransferEntropy); library(openxlsx) })

raw_results <- foreach(task = tasks, .combine = c, .errorhandling = "pass",
                       .packages = c("dplyr", "RTransferEntropy", "openxlsx")) %dopar% {
                         
                         res <- tryCatch({
                           process_single_group(
                             band1 = task$band1, band2 = task$band2, group = task$group,
                             subject_range = task$range, base_path = base_path, condition = condition,
                             microstates = microstates, lx = lx, ly = ly, nboot = nboot,
                             combo_output_dir = task$dir, combo_name = task$name
                           )
                         }, error = function(e) {
                           message("Error in task: ", task$name, " ", task$group)
                           return(NULL)
                         })
                         
                         rtn_list <- list(list(
                           combo_name = task$name,
                           group = task$group,
                           data = res
                         ))
                         return(rtn_list)
                       }

stopCluster(cl)
registerDoSEQ()
gc()

cat("\nParallel computation complete. Aggregating results...\n")

combo_map <- new.env(hash = TRUE)
for (combo in band_combinations) {
  cname <- paste0(combo[1], "-", combo[2])
  combo_map[[cname]] <- list(HC = list(), MDD = list())
}

for (item in raw_results) {
  if (is.null(item$data)) next
  cname <- item$combo_name
  grp <- item$group
  if (exists(cname, envir = combo_map)) {
    combo_map[[cname]][[grp]] <- item$data
  }
}

overall_stats <- list()
for (combo in band_combinations) {
  cname <- paste0(combo[1], "-", combo[2])
  combo_dir <- file.path(output_dir, cname)
  
  if (!exists(cname, envir = combo_map)) next
  
  hc_dat <- combo_map[[cname]]$HC
  mdd_dat <- combo_map[[cname]]$MDD
  
  cat("Generating report:", cname, "\n")
  
  stats <- finalize_band_combo(
    band1 = combo[1], band2 = combo[2],
    hc_results = hc_dat, mdd_results = mdd_dat,
    combo_output_dir = combo_dir, combo_name = cname
  )
  
  if (!is.null(stats)) {
    overall_stats[[length(overall_stats) + 1]] <- stats
  }
}

if (length(overall_stats) > 0) {
  overall_stats <- do.call(rbind, overall_stats)
  write.csv(overall_stats, file.path(output_dir, "overall_stats.csv"), row.names = F)
}

cat("\nGenerating comprehensive report...\n")
wb <- createWorkbook()
has_worksheet <- FALSE
if (length(overall_stats) > 0) {
  addWorksheet(wb, "Overall_Stats")
  writeData(wb, "Overall_Stats", overall_stats)
  has_worksheet <- TRUE
}

all_fdr <- list()
all_group_comp <- list()
for (combo in band_combinations) {
  name <- paste0(combo[1], "-", combo[2])
  dir <- file.path(output_dir, name)
  if (dir.exists(dir)) {
    f <- list.files(dir, paste0(name, "_FDR.csv"), full.names = T)
    if (length(f) > 0) { all_fdr[[name]] <- read.csv(f[1]); all_fdr[[name]]$band_combo <- name }
    f <- list.files(dir, paste0(name, "_group_diff.csv"), full.names = T)
    if (length(f) > 0) all_group_comp[[name]] <- read.csv(f[1])
  }
}

if (length(all_fdr) > 0) {
  addWorksheet(wb, "All_FDR")
  writeData(wb, "All_FDR", do.call(rbind, all_fdr))
  has_worksheet <- TRUE
}
if (length(all_group_comp) > 0) {
  addWorksheet(wb, "All_Group_Diff")
  writeData(wb, "All_Group_Diff", do.call(rbind, all_group_comp))
  has_worksheet <- TRUE
}
if (!has_worksheet) {
  addWorksheet(wb, "Info")
  writeData(wb, "Info", data.frame(Message = "No valid results"))
}

saveWorkbook(wb, file.path(output_dir, "transfer_entropy_report.xlsx"), overwrite = T)

end_time <- Sys.time()
cat("\n======== ANALYSIS COMPLETE ========\n")
cat("Total time:", round(difftime(end_time, start_time, units = "hours"), 2), "hours\n")