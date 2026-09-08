# ============================================================================
# Running the robust_shape_regression.R pipeline on the real shapes::rats data
#
# Your data layout (a "flat" format, NOT the usual 4-D array):
#   dim(rats$x) = 8 x 2 x 144   -> 144 = 18 rats x 8 measurement times
#   rats$no     : length 144, the rat ID for each slice along the 3rd dim
#   rats$time   : length 144, the time (day) for each slice along the 3rd dim
#
# This means array_to_configs() from the original script (which assumes an
# 8x2x8x18 array) does NOT apply directly here. Instead we use rats$no and
# rats$time to pull out the slices belonging to each rat and sort them by
# time.
# ============================================================================

library(shapes)
data(rats)

# Needed to write the results to an .xlsx file (install once if missing):
#   install.packages("openxlsx")
library(openxlsx)

# The file you uploaded earlier, containing weighted_procrustes,
# weighted_gpa, fit_shape_model, fit_real_data,
# run_real_data_contamination_study, etc.:
source("robust_shape_regression.R")

# ---------------------------------------------------------------------------
# 1) Extract the configurations for one specific rat from the 8 x 2 x 144 array
# ---------------------------------------------------------------------------
extract_rat <- function(rat_id) {
  idx <- which(rats$no == rat_id)
  if (length(idx) == 0) stop("rat_id not found: ", rat_id)
  ord <- order(rats$time[idx])        # make sure time is ascending
  idx <- idx[ord]

  X_list <- lapply(idx, function(i) rats$x[, , i])   # each element: 8 x 2
  times  <- rats$time[idx]

  list(times = times, X_list = X_list, rat_id = rat_id)
}

rat_ids <- sort(unique(rats$no))
cat("Number of rats:", length(rat_ids), " | IDs:", paste(rat_ids, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# 2) Example: fit on one rat (e.g. the first rat ID present in the data)
# ---------------------------------------------------------------------------
r1 <- extract_rat(rat_ids[1])
str(r1$X_list[[1]])   # should be 8 x 2
print(r1$times)       # should be 7 14 21 30 40 60 90 150

fit1 <- fit_real_data(r1$times, r1$X_list, method = "robust")
print(fit1)

# ---------------------------------------------------------------------------
# 3) Contamination study on that same rat
#    (naive / aware / robust, matching Table 1 of the paper)
# ---------------------------------------------------------------------------
study1 <- run_real_data_contamination_study(r1$times, r1$X_list, n_rep = 40)
print(study1, row.names = FALSE)

# ---------------------------------------------------------------------------
# 4) Run the full pipeline over all 18 rats and aggregate results
#    (use a smaller n_rep, e.g. 10, first to check everything runs)
# ---------------------------------------------------------------------------
run_all_rats <- function(n_rep = 40) {
  all_results <- vector("list", length(rat_ids))
  names(all_results) <- as.character(rat_ids)

  for (id in rat_ids) {
    cat("== Rat", id, "==\n")
    d <- extract_rat(id)
    res <- run_real_data_contamination_study(d$times, d$X_list, n_rep = n_rep)
    res$RatID <- id
    all_results[[as.character(id)]] <- res
  }

  do.call(rbind, all_results)
}

# Full run, n_rep = 40 (matches the replication count used for Table 3 of
# the paper; takes a few minutes for all 18 rats). The updated
# run_real_data_contamination_study() now also returns, per row: a median
# and 10%-trimmed-mean version of R2/RMSD/Sim (robust to the strong right
# skew R2 shows under contamination -- see referee report, Major point 3),
# plus paired-Wilcoxon p-values (p_R2_vs_naive, p_R2_vs_aware) computed
# from methods evaluated on the SAME contaminated draw within each
# replicate (referee report, Major point 5). Use this table, not a
# n_rep = 5 "quick" run, for anything reported in the paper: the earlier
# n_rep = 5 pass that produced rats_results.xlsx was explicitly flagged as
# too noisy to be conclusive (see Section 5.3 / Appendix B of the paper).
final_table <- run_all_rats(n_rep = 40)
print(final_table, row.names = FALSE)
# write.csv(final_table, "rats_contamination_results.csv", row.names = FALSE)

# A much faster n_rep = 5 pass is still useful as a smoke test that the
# pipeline runs end to end on all 18 rats, but should NOT be used for any
# number quoted in the paper (its results are noticeably noisier -- e.g.
# it can flip the sign of small naive-vs-robust differences).
# final_table_smoketest <- run_all_rats(n_rep = 5)

# ---------------------------------------------------------------------------
# 4b) Cross-rat heterogeneity summary (Appendix B of the paper): aggregate
# the per-rat table above into one row per scenario/method, reporting the
# across-rat mean/min/max of Sim and R2 (scale-free, so they aggregate
# directly), and RMSD expressed as a ratio to each rat's own Clean/naive
# RMSD (since raw RMSD is in whatever digitization units that rat's
# configurations use, and is not directly comparable across rats -- see
# referee report / paper Section on heterogeneity for the full rationale).
# ---------------------------------------------------------------------------
summarize_rat_heterogeneity <- function(final_table) {
  clean_naive_rmsd <- setNames(
    final_table$RMSD_mean[final_table$Scenario == "Clean" & final_table$Method == "naive"],
    final_table$RatID[final_table$Scenario == "Clean" & final_table$Method == "naive"]
  )
  final_table$RMSD_ratio <- final_table$RMSD_mean / clean_naive_rmsd[as.character(final_table$RatID)]

  agg <- aggregate(
    cbind(Sim_mean, R2_mean, RMSD_ratio) ~ Scenario + Method,
    data = final_table,
    FUN = function(x) c(mean = mean(x), min = min(x), max = max(x))
  )
  agg
}

heterogeneity_table <- summarize_rat_heterogeneity(final_table)
print(heterogeneity_table, row.names = FALSE)

# ---------------------------------------------------------------------------
# 5) Save every result table to a single Excel workbook, one sheet each
# ---------------------------------------------------------------------------
save_results_to_excel <- function(fit1, study1, final_table,
                                   file = "rats_results.xlsx") {
  wb <- createWorkbook()

  # Sheet 1: robust fit on rat_ids[1] (Xtilde is a list of matrices, so we
  # stack it into one long data frame with a Time / Landmark column)
  addWorksheet(wb, "Fit_Rat1")
  xt <- fit1$Xtilde
  fit_df <- do.call(rbind, lapply(seq_along(xt), function(i) {
    m <- as.data.frame(xt[[i]])
    colnames(m) <- paste0("Coord", seq_len(ncol(m)))
    cbind(Time = r1$times[i], Landmark = seq_len(nrow(m)), m)
  }))
  writeData(wb, "Fit_Rat1", fit_df)

  # Sheet 2: contamination study on rat_ids[1] (naive/aware/robust x scenario)
  addWorksheet(wb, "Contamination_Rat1")
  writeData(wb, "Contamination_Rat1", study1)

  # Sheet 3: contamination study for every rat individually (n_rep = 40),
  # now including median/trimmed-mean columns and paired-test p-values
  addWorksheet(wb, "Contamination_AllRats")
  writeData(wb, "Contamination_AllRats", final_table)

  # Sheet 4: cross-rat heterogeneity summary (Appendix B of the paper)
  addWorksheet(wb, "Heterogeneity_Summary")
  writeData(wb, "Heterogeneity_Summary", heterogeneity_table)

  saveWorkbook(wb, file, overwrite = TRUE)
  cat("Saved:", normalizePath(file), "\n")
}

save_results_to_excel(fit1, study1, final_table, file = "rats_results.xlsx")
