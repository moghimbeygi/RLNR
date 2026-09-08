# ============================================================================
# Re-run EVERYTHING with the updated robust_shape_regression.R and n_rep=40:
#   1) Table 2  -- hexagon simulation study
#   2) Table 3  -- single representative rat, real-data contamination study
#   3) Appendix B -- all 18 rats, real-data contamination study
#
# Run this on a machine with R + the "shapes" package installed. This
# environment (the one used to prepare the paper/code) has neither R nor
# internet access, so none of the numbers below have actually been
# generated yet -- this script is what produces them.
#
# Usage:
#   Rscript rerun_all_n40.R
# or from an interactive R session:
#   source("rerun_all_n40.R")
#
# Everything is saved to rerun_results.xlsx (one sheet per table) and to
# individual .csv files, so you can either hand the .xlsx back for the
# tables to be transcribed into the paper, or paste the console output.
# ============================================================================

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx", repos = "https://cloud.r-project.org")
if (!requireNamespace("shapes",   quietly = TRUE)) install.packages("shapes",   repos = "https://cloud.r-project.org")

library(shapes)
library(openxlsx)
data("rats")

# Path to the UPDATED script (the one with contamination_scale(),
# robust_summary(), paired_method_test(), and the common-random-numbers
# restructuring of run_simulation()/run_real_data_contamination_study()).
# Adjust this path if the file lives elsewhere on your machine.
source("robust_shape_regression.R")

set.seed(1)

# ---------------------------------------------------------------------------
# 1) Table 2: hexagon simulation, n_rep = 40
# ---------------------------------------------------------------------------
cat("\n============ Table 2: hexagon simulation (n_rep = 40) ============\n")
table2 <- run_simulation(n_rep = 40)
print(table2, row.names = FALSE)

# ---------------------------------------------------------------------------
# 2) Table 3: single representative rat, n_rep = 40
#    (uses the same rat used originally for Table 3; adjust rat_index if
#    your original Table 3 used a different one -- the classic illustrative
#    choice in most shapes::rats analyses is the 2nd rat)
# ---------------------------------------------------------------------------
extract_rat <- function(rat_id) {
  idx <- which(rats$no == rat_id)
  if (length(idx) == 0) stop("rat_id not found: ", rat_id)
  ord <- order(rats$time[idx])
  idx <- idx[ord]
  list(times = rats$time[idx],
       X_list = lapply(idx, function(i) rats$x[, , i]),
       rat_id = rat_id)
}

rat_ids <- sort(unique(rats$no))
rat_index <- 2                      # <-- change this if Table 3 used a different rat
r_table3 <- extract_rat(rat_ids[rat_index])

cat("\n============ Table 3: single rat (ID =", r_table3$rat_id, "), n_rep = 40 ============\n")
table3 <- run_real_data_contamination_study(r_table3$times, r_table3$X_list, n_rep = 40)
print(table3, row.names = FALSE)

# ---------------------------------------------------------------------------
# 3) Appendix B: all 18 rats, n_rep = 40 each
# ---------------------------------------------------------------------------
cat("\n============ Appendix B: all 18 rats, n_rep = 40 each ============\n")
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
appendixB_allrats <- run_all_rats(n_rep = 40)

# Cross-rat heterogeneity summary (Sim/R2 aggregate directly; RMSD is
# expressed as a ratio to each rat's own Clean/naive RMSD, since raw RMSD
# units need not agree across rats -- see Appendix B of the paper)
summarize_rat_heterogeneity <- function(final_table) {
  clean_naive_rmsd <- setNames(
    final_table$RMSD_mean[final_table$Scenario == "Clean" & final_table$Method == "naive"],
    final_table$RatID[final_table$Scenario == "Clean" & final_table$Method == "naive"]
  )
  final_table$RMSD_ratio <- final_table$RMSD_mean / clean_naive_rmsd[as.character(final_table$RatID)]
  aggregate(
    cbind(Sim_mean, R2_mean, RMSD_ratio) ~ Scenario + Method,
    data = final_table,
    FUN = function(x) c(mean = mean(x), min = min(x), max = max(x))
  )
}
appendixB_heterogeneity <- summarize_rat_heterogeneity(appendixB_allrats)
print(appendixB_heterogeneity, row.names = FALSE)

# ---------------------------------------------------------------------------
# 4) Save everything
# ---------------------------------------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "Table2_Hexagon");            writeData(wb, "Table2_Hexagon", table2)
addWorksheet(wb, "Table3_SingleRat");          writeData(wb, "Table3_SingleRat", table3)
addWorksheet(wb, "AppendixB_AllRats");         writeData(wb, "AppendixB_AllRats", appendixB_allrats)
addWorksheet(wb, "AppendixB_Heterogeneity");   writeData(wb, "AppendixB_Heterogeneity", appendixB_heterogeneity)
saveWorkbook(wb, "rerun_results.xlsx", overwrite = TRUE)

write.csv(table2, "table2_hexagon_n40.csv", row.names = FALSE)
write.csv(table3, "table3_singlerat_n40.csv", row.names = FALSE)
write.csv(appendixB_allrats, "appendixB_allrats_n40.csv", row.names = FALSE)
write.csv(appendixB_heterogeneity, "appendixB_heterogeneity_n40.csv", row.names = FALSE)

cat("\nSaved: rerun_results.xlsx (4 sheets) and 4 individual .csv files.\n")
cat("Send rerun_results.xlsx back and the paper's Tables 2, 3, and Appendix B\n")
cat("will be updated with these numbers (including the new median/trimmed-mean\n")
cat("columns and paired-test p-values).\n")

