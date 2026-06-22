if (!require(dplyr)) install.packages("dplyr")
if (!require(reshape2)) install.packages("reshape2")
library(dplyr)
library(reshape2)

base_path <- "path/to/EC_sequences"
output_dir <- "path/to/output/cooccurrence"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

hc_range <- 1:162
mdd_range <- 6:416
condition <- "EC"
microstates <- c("A", "B", "C", "D")

band_combinations <- list(
  c("Broadband", "Delta"),
  c("Broadband", "Theta"), 
  c("Broadband", "Alpha"),
  c("Broadband", "Beta"),
  c("Broadband", "Gamma")
)

process_single_subject_bands <- function(subject_id, group, base_path, condition, band1, band2) {
  
  band1_lower <- tolower(band1)
  band2_lower <- tolower(band2)
  
  band1_file <- file.path(base_path, band1_lower, 
                          paste0(group, subject_id, "_", condition, "_", band1_lower, ".set_dynamics.txt"))
  band2_file <- file.path(base_path, band2_lower,
                          paste0(group, subject_id, "_", condition, "_", band2_lower, ".set_dynamics.txt"))
  
  if (!file.exists(band1_file) | !file.exists(band2_file)) {
    return(NULL)
  }
  
  data_band1 <- read.table(
    file = band1_file,
    header = TRUE,
    sep = "",
    strip.white = TRUE,
    na.strings = c("NA", ""),
    stringsAsFactors = FALSE,
    encoding = "UTF-8"
  )
  
  data_band2 <- read.table(
    file = band2_file,
    header = TRUE,
    sep = "",
    strip.white = TRUE,
    na.strings = c("NA", ""),
    stringsAsFactors = FALSE,
    encoding = "UTF-8"
  )
  
  convert_to_labels <- function(data) {
    intensity_cols <- c("A", "B", "C", "D")
    
    labels <- apply(data[, intensity_cols], 1, function(row) {
      abs_values <- abs(row)
      max_abs <- max(abs_values)
      
      if (max_abs == 0) {
        return(NA)
      } else {
        return(intensity_cols[which.max(abs_values)])
      }
    })
    
    return(data.frame(
      Time = data$Time,
      Label = labels,
      stringsAsFactors = FALSE
    ))
  }
  
  labels_band1 <- convert_to_labels(data_band1)
  labels_band2 <- convert_to_labels(data_band2)
  
  merged_data <- merge(labels_band1, labels_band2, 
                       by = "Time", 
                       suffixes = c(paste0("_", band1), paste0("_", band2)))
  colnames(merged_data) <- c("Time", band1, band2)
  
  valid_data <- na.omit(merged_data)
  
  if (nrow(valid_data) == 0) {
    return(NULL)
  }
  
  cooccurrence_freq <- table(valid_data[[band1]], valid_data[[band2]])
  
  cooccurrence_prob_row <- prop.table(cooccurrence_freq, margin = 1)
  cooccurrence_prob_col <- prop.table(cooccurrence_freq, margin = 2)
  cooccurrence_prob_joint <- prop.table(cooccurrence_freq)
  
  marginal_prob_band1 <- prop.table(table(valid_data[[band1]]))
  marginal_prob_band2 <- prop.table(table(valid_data[[band2]]))
  
  deltaP_prob_row <- matrix(NA, nrow = 4, ncol = 4)
  rownames(deltaP_prob_row) <- microstates
  colnames(deltaP_prob_row) <- microstates
  
  deltaP_prob_col <- matrix(NA, nrow = 4, ncol = 4)
  rownames(deltaP_prob_col) <- microstates
  colnames(deltaP_prob_col) <- microstates
  
  for (b1_state in microstates) {
    for (b2_state in microstates) {
      if (b1_state %in% rownames(cooccurrence_prob_row) && 
          b2_state %in% colnames(cooccurrence_prob_row)) {
        raw_prob_row <- cooccurrence_prob_row[b1_state, b2_state]
        baseline_row <- ifelse(b2_state %in% names(marginal_prob_band2), 
                               marginal_prob_band2[b2_state], 0)
        deltaP_prob_row[b1_state, b2_state] <- raw_prob_row - baseline_row
      }
      
      if (b1_state %in% rownames(cooccurrence_prob_col) && 
          b2_state %in% colnames(cooccurrence_prob_col)) {
        raw_prob_col <- cooccurrence_prob_col[b1_state, b2_state]
        baseline_col <- ifelse(b1_state %in% names(marginal_prob_band1), 
                               marginal_prob_band1[b1_state], 0)
        deltaP_prob_col[b1_state, b2_state] <- raw_prob_col - baseline_col
      }
    }
  }
  
  extract_all_pairs <- function(prob_matrix, prefix) {
    pairs <- c()
    for (b1_state in microstates) {
      for (b2_state in microstates) {
        prob_value <- ifelse(b1_state %in% rownames(prob_matrix) && 
                               b2_state %in% colnames(prob_matrix),
                             prob_matrix[b1_state, b2_state], NA)
        pair_name <- paste0(prefix, "_", b1_state, b2_state)
        pairs[pair_name] <- prob_value
      }
    }
    return(pairs)
  }
  
  row_pairs <- extract_all_pairs(cooccurrence_prob_row, "row")
  col_pairs <- extract_all_pairs(cooccurrence_prob_col, "col")
  joint_pairs <- extract_all_pairs(cooccurrence_prob_joint, "joint")
  deltaP_row_pairs <- extract_all_pairs(deltaP_prob_row, "deltaP_row")
  deltaP_col_pairs <- extract_all_pairs(deltaP_prob_col, "deltaP_col")
  
  marginal_band1 <- c()
  marginal_band2 <- c()
  for (state in microstates) {
    marginal_band1[paste0("marginal_", band1, "_", state)] <- 
      ifelse(state %in% names(marginal_prob_band1), marginal_prob_band1[state], NA)
    marginal_band2[paste0("marginal_", band2, "_", state)] <- 
      ifelse(state %in% names(marginal_prob_band2), marginal_prob_band2[state], NA)
  }
  
  all_pairs <- c(row_pairs, col_pairs, joint_pairs, 
                 deltaP_row_pairs, deltaP_col_pairs,
                 marginal_band1, marginal_band2)
  
  result <- list(
    subject_id = paste0(group, subject_id),
    group = group,
    condition = condition,
    band1 = band1,
    band2 = band2,
    n_valid = nrow(valid_data),
    all_probabilities = all_pairs
  )
  
  return(result)
}

for (band_combo in band_combinations) {
  band1 <- band_combo[1]
  band2 <- band_combo[2]
  combo_name <- paste0(band1, "-", band2)
  
  combo_output_dir <- file.path(output_dir, combo_name)
  if (!dir.exists(combo_output_dir)) {
    dir.create(combo_output_dir, recursive = TRUE)
  }
  
  combo_results <- list()
  
  hc_count <- 0
  for (i in hc_range) {
    result <- process_single_subject_bands(i, "HC", base_path, condition, band1, band2)
    if (!is.null(result)) {
      combo_results[[length(combo_results) + 1]] <- result
      hc_count <- hc_count + 1
    }
  }
  
  mdd_count <- 0
  for (i in mdd_range) {
    result <- process_single_subject_bands(i, "MDD", base_path, condition, band1, band2)
    if (!is.null(result)) {
      combo_results[[length(combo_results) + 1]] <- result
      mdd_count <- mdd_count + 1
    }
  }
  
  if (length(combo_results) == 0) {
    next
  }
  
  detailed_df <- do.call(rbind, lapply(combo_results, function(x) {
    base_info <- data.frame(
      subject_id = x$subject_id,
      group = x$group,
      condition = x$condition,
      band1 = x$band1,
      band2 = x$band2,
      n_valid = x$n_valid,
      stringsAsFactors = FALSE
    )
    
    prob_df <- as.data.frame(t(x$all_probabilities))
    cbind(base_info, prob_df)
  }))
  
  rownames(detailed_df) <- NULL
  
  probability_types <- c("row", "col", "joint", "deltaP_row", "deltaP_col")
  
  summary_tables <- list()
  
  for (prob_type in probability_types) {
    prob_cols <- grep(paste0("^", prob_type, "_"), names(detailed_df), value = TRUE)
    
    if (length(prob_cols) > 0) {
      summary_table <- detailed_df %>%
        group_by(group) %>%
        summarise(across(all_of(prob_cols), 
                         list(mean = ~mean(., na.rm = TRUE),
                              sd = ~sd(., na.rm = TRUE)),
                         .names = "{.col}_{.fn}"))
      
      summary_tables[[prob_type]] <- summary_table
    }
  }
  
  marginal_cols_band1 <- grep(paste0("^marginal_", band1), names(detailed_df), value = TRUE)
  marginal_cols_band2 <- grep(paste0("^marginal_", band2), names(detailed_df), value = TRUE)
  
  if (length(marginal_cols_band1) > 0) {
    marginal_summary_band1 <- detailed_df %>%
      group_by(group) %>%
      summarise(across(all_of(marginal_cols_band1), 
                       list(mean = ~mean(., na.rm = TRUE),
                            sd = ~sd(., na.rm = TRUE)),
                       .names = "{.col}_{.fn}"))
    summary_tables[["marginal_band1"]] <- marginal_summary_band1
  }
  
  if (length(marginal_cols_band2) > 0) {
    marginal_summary_band2 <- detailed_df %>%
      group_by(group) %>%
      summarise(across(all_of(marginal_cols_band2), 
                       list(mean = ~mean(., na.rm = TRUE),
                            sd = ~sd(., na.rm = TRUE)),
                       .names = "{.col}_{.fn}"))
    summary_tables[["marginal_band2"]] <- marginal_summary_band2
  }
  
  write.csv(detailed_df, file.path(combo_output_dir, paste0(combo_name, "_microstate_cooccurrence_detailed.csv")), 
            row.names = FALSE, fileEncoding = "UTF-8")
  
  for (prob_type in names(summary_tables)) {
    if (!is.null(summary_tables[[prob_type]])) {
      write.csv(summary_tables[[prob_type]], 
                file.path(combo_output_dir, paste0(combo_name, "_", prob_type, "_summary.csv")), 
                row.names = FALSE, fileEncoding = "UTF-8")
    }
  }
  
  saveRDS(combo_results, file.path(combo_output_dir, paste0(combo_name, "_cooccurrence_results.rds")))
}