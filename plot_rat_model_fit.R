# ============================================================================
# Extends plot_rat_growth.R: overlays the paper's model (naive/aware/robust)
# on top of the raw observed skull landmarks for one chosen rat, in two ways:
#
#   1) plot_rat_fit_overlay()      -- per age (faceted), the OBSERVED
#                                      configuration vs. the model's
#                                      leave-one-out PREDICTED configuration,
#                                      drawn on the same axes
#   2) plot_rat_pred_vs_actual()   -- a predicted-vs-actual scatter over every
#                                      (age, landmark, coordinate), the same
#                                      comparison that produces the R2/RMSD/Sim
#                                      numbers in the paper's tables, so you
#                                      can see directly what those numbers
#                                      are summarizing
#
# Requires: shapes, ggplot2, and robust_shape_regression.R in the same folder.
# ============================================================================

library(shapes)
library(ggplot2)

data(rats)
source("robust_shape_regression.R")   # fit_shape_model(), nw_predict_shape(), eval_scenario()

# ---------------------------------------------------------------------------
# Same helpers as plot_rat_growth.R (repeated here so this file is
# self-contained)
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

rat_to_df <- function(rat_id) {
  d <- extract_rat(rat_id)
  k <- nrow(d$X_list[[1]])
  do.call(rbind, lapply(seq_along(d$times), function(i) {
    data.frame(RatID = rat_id, Time = d$times[i], Landmark = seq_len(k),
               X = d$X_list[[i]][, 1], Y = d$X_list[[i]][, 2])
  }))
}

# ---------------------------------------------------------------------------
# Fit the model (naive / aware / robust) and compute the SAME leave-one-out
# predictions used by the paper's evaluation. This mirrors the bandwidth
# grid-search inside run_method_once() in robust_shape_regression.R exactly,
# but returns the fitted alignment (Xtilde) and the LOO predictions
# themselves, rather than only the summary R2/RMSD/Sim numbers.
# ---------------------------------------------------------------------------
fit_and_predict_loo <- function(times, X_list, method = "robust", h_grid = NULL) {
  n <- length(times); k <- nrow(X_list[[1]])
  delta_list <- replicate(n, rep(1, k), simplify = FALSE)   # real data here: nothing missing

  if (is.null(h_grid)) {
    span <- diff(range(times))
    h_grid <- seq(span / 40, span / 2, length.out = 15)
  }

  fitted <- fit_shape_model(times, X_list, delta_list, robust = (method == "robust"))
  Xtilde <- fitted$Xtilde
  W      <- fitted$weights

  best_h <- h_grid[1]; best_err <- Inf
  for (h in h_grid) {
    errs <- numeric(0)
    for (i in seq_len(n)) {
      p_i <- nw_predict_shape(times[i], times, Xtilde, W, h, exclude = i)
      wl  <- W[[i]]
      if (!anyNA(p_i) && sum(wl) > 1e-8) {
        errs <- c(errs, sum(wl * rowSums((p_i - Xtilde[[i]])^2)) / sum(wl))
      }
    }
    if (length(errs) > 0) {
      err <- stats::median(errs)
      if (err < best_err) { best_err <- err; best_h <- h }
    }
  }

  preds <- lapply(seq_len(n), function(i)
    nw_predict_shape(times[i], times, Xtilde, W, best_h, exclude = i))

  mean_true <- Reduce(`+`, X_list) / n
  metrics <- eval_scenario(preds, X_list, mean_true)   # same R2/RMSD/Sim as the paper's tables

  list(Xtilde = Xtilde, preds = preds, best_h = best_h, metrics = metrics)
}

# ---------------------------------------------------------------------------
# 1) Per-age overlay: observed configuration vs. leave-one-out prediction
# ---------------------------------------------------------------------------
plot_rat_fit_overlay <- function(rat_id, method = "robust") {
  d <- extract_rat(rat_id)
  fit <- fit_and_predict_loo(d$times, d$X_list, method = method)
  k <- nrow(d$X_list[[1]])

  df_obs <- do.call(rbind, lapply(seq_along(d$times), function(i)
    data.frame(Time = d$times[i], Landmark = seq_len(k),
               X = d$X_list[[i]][, 1], Y = d$X_list[[i]][, 2], Source = "Observed")))

  df_pred <- do.call(rbind, lapply(seq_along(d$times), function(i) {
    p <- fit$preds[[i]]
    if (anyNA(p)) return(NULL)   # LOO prediction undefined at that age/bandwidth
    data.frame(Time = d$times[i], Landmark = seq_len(k),
               X = p[, 1], Y = p[, 2], Source = "LOO prediction")
  }))

  df <- rbind(df_obs, df_pred)
  df$Time   <- factor(df$Time, levels = sort(unique(d$times)),
                       labels = paste0(sort(unique(d$times)), "d"))
  df$Source <- factor(df$Source, levels = c("Observed", "LOO prediction"))

  ggplot(df, aes(x = X, y = Y, group = interaction(Time, Source),
                  color = Source, linetype = Source, fill = Source)) +
    geom_polygon(alpha = 0.12, linewidth = 0.8) +
    geom_point(size = 1.8) +
    scale_color_manual(values = c("Observed" = "steelblue4", "LOO prediction" = "firebrick")) +
    scale_fill_manual(values  = c("Observed" = "steelblue4", "LOO prediction" = "firebrick")) +
    scale_linetype_manual(values = c("Observed" = "solid", "LOO prediction" = "dashed")) +
    coord_fixed() +
    facet_wrap(~ Time, nrow = 2) +
    labs(
      title    = paste0("Rat ", rat_id, ": observed vs. ", method, " leave-one-out prediction"),
      subtitle = sprintf("R2 = %.1f, RMSD = %.2f, Sim = %.1f (h = %.3f)",
                          fit$metrics$R2, fit$metrics$RMSD, fit$metrics$Sim, fit$best_h),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          strip.text = element_text(face = "bold"), legend.position = "bottom")
}

# ---------------------------------------------------------------------------
# 2) Predicted-vs-actual scatter: every (age, landmark, coordinate) triple,
#    Procrustes-aligned prediction vs. the true observed value -- literally
#    the comparison that R2/RMSD/Sim in the paper's tables summarize.
# ---------------------------------------------------------------------------
plot_rat_pred_vs_actual <- function(rat_id, method = "robust") {
  d <- extract_rat(rat_id)
  fit <- fit_and_predict_loo(d$times, d$X_list, method = method)
  k <- nrow(d$X_list[[1]])

  rows <- list()
  for (i in seq_along(d$times)) {
    p <- fit$preds[[i]]
    if (anyNA(p)) next
    # align the prediction onto the true configuration before comparing,
    # exactly as eval_scenario() does internally, so the scatter reflects
    # shape error rather than an arbitrary rotation/translation mismatch
    aligned <- procrustes_align_full(p, d$X_list[[i]])
    rows[[length(rows) + 1]] <- data.frame(
      Time = d$times[i], Landmark = seq_len(k),
      Actual_X = d$X_list[[i]][, 1], Pred_X = aligned[, 1],
      Actual_Y = d$X_list[[i]][, 2], Pred_Y = aligned[, 2]
    )
  }
  df <- do.call(rbind, rows)
  df_long <- rbind(
    data.frame(Time = df$Time, Landmark = df$Landmark, Coord = "X",
               Actual = df$Actual_X, Predicted = df$Pred_X),
    data.frame(Time = df$Time, Landmark = df$Landmark, Coord = "Y",
               Actual = df$Actual_Y, Predicted = df$Pred_Y)
  )

  rng <- range(c(df_long$Actual, df_long$Predicted))

  ggplot(df_long, aes(x = Actual, y = Predicted, color = factor(Time))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(size = 2, alpha = 0.85) +
    coord_fixed(xlim = rng, ylim = rng) +
    facet_wrap(~ Coord, labeller = label_both) +
    scale_color_viridis_d(name = "Age (days)") +
    labs(
      title    = paste0("Rat ", rat_id, ": ", method, " leave-one-out prediction vs. actual"),
      subtitle = sprintf("R2 = %.1f, RMSD = %.2f, Sim = %.1f  (dashed line = perfect prediction)",
                          fit$metrics$R2, fit$metrics$RMSD, fit$metrics$Sim),
      x = "Actual coordinate", y = "Predicted coordinate (aligned)"
    ) +
    theme_minimal(base_size = 13)
}

# ---------------------------------------------------------------------------
# Run it
# ---------------------------------------------------------------------------
rat_id <- 2
method <- "robust"   # "naive", "aware", or "robust"

p_overlay_fit <- plot_rat_fit_overlay(rat_id, method)
p_pred_actual <- plot_rat_pred_vs_actual(rat_id, method)

print(p_overlay_fit)
print(p_pred_actual)

ggsave(sprintf("rat%d_%s_fit_overlay.png",  rat_id, method), p_overlay_fit, width = 9, height = 5,   dpi = 150)
ggsave(sprintf("rat%d_%s_pred_vs_actual.png", rat_id, method), p_pred_actual, width = 8, height = 4.5, dpi = 150)
