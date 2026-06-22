library(infotheo)
library(aricode)

base_path <- "path/to/EC_sequences"
output_dir <- "path/to/output/AMI"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

bands <- c("Broadband", "Delta", "Theta", "Alpha", "Beta", "Gamma")
hc_ids <- 1:162
mdd_ids <- 6:416
all_ids <- list(HC = hc_ids, MDD = mdd_ids)

read_aligned_labels <- function(group, subject_id, bands, condition = "EC") {
  raw_list <- list()
  for (b in bands) {
    band_lower <- tolower(b)
    file_path <- file.path(base_path, band_lower,
                           paste0(group, subject_id, "_", condition, "_", band_lower, ".set_dynamics.txt"))
    if (!file.exists(file_path)) return(NULL)
    raw_list[[b]] <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE)
  }
  
  time_sets <- lapply(raw_list, function(df) df$Time)
  common_times <- Reduce(intersect, time_sets)
  if (length(common_times) == 0) return(NULL)
  
  label_mat <- matrix(NA_character_, nrow = length(common_times), ncol = length(bands))
  colnames(label_mat) <- bands
  for (b in bands) {
    df <- raw_list[[b]]
    idx <- match(common_times, df$Time)
    sub <- df[idx, c("A", "B", "C", "D")]
    labs <- apply(sub, 1, function(row) {
      abs_vals <- abs(row)
      if (max(abs_vals) == 0) return(NA_character_)
      c("A", "B", "C", "D")[which.max(abs_vals)]
    })
    label_mat[, b] <- labs
  }
  
  valid <- complete.cases(label_mat)
  if (sum(valid) == 0) return(NULL)
  
  result <- list()
  for (b in bands) {
    result[[b]] <- factor(label_mat[valid, b], levels = c("A", "B", "C", "D"))
  }
  return(result)
}

label_seqs <- list()
for (b in bands) label_seqs[[b]] <- list()

for (grp in c("HC", "MDD")) {
  ids <- all_ids[[grp]]
  for (sid in ids) {
    aligned <- read_aligned_labels(grp, sid, bands, "EC")
    if (is.null(aligned)) next
    subj_key <- paste0(grp, sid)
    for (b in bands) {
      label_seqs[[b]][[subj_key]] <- aligned[[b]]
    }
  }
}

compute_ami <- function(x, y) {
  AMI(x, y)
}

common_subjects <- Reduce(intersect, lapply(label_seqs, names))

K <- length(bands)
ami_matrix <- matrix(NA, K, K, dimnames = list(bands, bands))

ami_cache <- list()
cache_idx <- 1

for (i in 1:K) {
  for (j in i:K) {
    band_i <- bands[i]
    band_j <- bands[j]
    if (i == j) {
      ami_matrix[i, j] <- 1
      next
    }
    ami_values <- sapply(common_subjects, function(subj) {
      xi <- label_seqs[[band_i]][[subj]]
      xj <- label_seqs[[band_j]][[subj]]
      compute_ami(xi, xj)
    })
    mean_ami <- mean(ami_values, na.rm = TRUE)
    ami_matrix[i, j] <- mean_ami
    ami_matrix[j, i] <- mean_ami
    
    ami_cache[[cache_idx]] <- list(b1 = band_i, b2 = band_j,
                                   ami_per_subj = ami_values)
    cache_idx <- cache_idx + 1
  }
}

write.csv(ami_matrix, file.path(output_dir, "AMI_matrix_EC.csv"), row.names = TRUE)

set.seed(123)
n_perm <- 1000
n_subj <- length(common_subjects)

p_raw_matrix <- matrix(NA, K, K, dimnames = list(bands, bands))
p_fdr_matrix <- matrix(NA, K, K, dimnames = list(bands, bands))

total_pairs <- K * (K - 1) / 2
pair_done <- 0
p_vals_upper <- c()
pair_indices <- list()

for (i in 1:(K - 1)) {
  for (j in (i + 1):K) {
    band_i <- bands[i]
    band_j <- bands[j]
    pair_done <- pair_done + 1
    
    cache_entry <- Filter(function(x) x$b1 == band_i && x$b2 == band_j,
                          ami_cache)[[1]]
    ami_obs <- cache_entry$ami_per_subj
    obs_mean <- mean(ami_obs)
    
    null_means <- replicate(n_perm, {
      null_ami <- sapply(common_subjects, function(subj) {
        xi <- label_seqs[[band_i]][[subj]]
        xj_shuffled <- sample(label_seqs[[band_j]][[subj]])
        compute_ami(xi, xj_shuffled)
      })
      mean(null_ami)
    })
    
    p_raw <- (sum(null_means >= obs_mean) + 1) / (n_perm + 1)
    
    p_raw_matrix[i, j] <- p_raw
    p_raw_matrix[j, i] <- p_raw
    p_vals_upper <- c(p_vals_upper, p_raw)
    pair_indices[[length(pair_indices) + 1]] <- c(i, j)
  }
}

p_fdr_upper <- p.adjust(p_vals_upper, method = "fdr")
for (idx in seq_along(pair_indices)) {
  i <- pair_indices[[idx]][1]
  j <- pair_indices[[idx]][2]
  p_fdr_matrix[i, j] <- p_fdr_upper[idx]
  p_fdr_matrix[j, i] <- p_fdr_upper[idx]
}

write.csv(p_raw_matrix, file.path(output_dir, "AMI_permutation_pvalues_raw.csv"), row.names = TRUE)
write.csv(p_fdr_matrix, file.path(output_dir, "AMI_permutation_pvalues_fdr.csv"), row.names = TRUE)