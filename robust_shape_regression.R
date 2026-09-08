# ============================================================================
# Robust, Missing-Data-Aware Nonparametric Longitudinal Regression for Shape
# Data (companion implementation)
#
# Implements:
#   - weighted (diagonal, landmark-specific) Procrustes alignment via SVD, eq (5)
#   - two-pass pilot -> Tukey biweight reweighting scheme, Section 3.2
#   - doubly weighted Nadaraya-Watson shape predictor, eq (7)
#   - hexagon growth simulation + contamination scenarios of Table 1
#   - naive / aware / robust estimators and Table-2-style evaluation
#
# Revisions made in response to referee comments on the accompanying
# manuscript (see the referee report shared alongside this file):
#   - explicit, documented outlier-magnitude formula for real data
#     (contamination_scale()); removed a previously dead/unused
#     intermediate calculation that did not match what was actually applied
#   - robust_summary(): median and 10%-trimmed mean reported alongside the
#     mean for every criterion, since R2 in particular can be strongly
#     right-skewed under contamination
#   - naive / aware / robust are now evaluated on a COMMON contaminated
#     draw within each replicate ("common random numbers"), rather than
#     each method independently re-randomizing; this reduces Monte Carlo
#     noise in the comparison and is required for the paired test below
#   - paired_method_test(): a paired Wilcoxon signed-rank test (robust vs.
#     naive, robust vs. aware) reported alongside every summary table
#   - removed stray non-R console output that had been appended to the end
#     of this file, which would break source() on this file as previously
#     distributed
#
# Further revisions made in response to a second referee report reviewing
# this file directly:
#   - fit_shape_model() / fit_shape_model_config_level(): the robust scale
#     estimate sigma_hat is now calibrated against the chi_p distribution
#     of ||e_{i,l}|| (p = ambient dimension), not the classical 1-D MAD
#     constant 1.4826, which assumes a scalar Gaussian deviation and does
#     not apply to a p-dimensional residual norm. For p=2 (all simulation
#     and real-data landmarks in this paper) the old formula overstated
#     sigma_hat by a factor of about 1.75, loosening the Tukey biweight
#     cutoff c*sigma_hat well past its nominal 95%-efficiency calibration.
#   - paired_method_test(): now also returns NA when both compared series
#     are individually constant across "replications" (e.g. a Clean
#     real-data scenario that draws no random contamination at all), since
#     a paired test on non-independent, non-random repeats is not a valid
#     inference even when the two series differ from each other.
#
# Base R only, no external packages required.
# ============================================================================

# ---------------------------------------------------------------------------
# 1. Weighted orthogonal Procrustes problem, closed-form SVD solution (eq 5)
# ---------------------------------------------------------------------------
#' @param X k x p source configuration
#' @param Y k x p target configuration
#' @param w length-k nonnegative weight vector (0 = missing / fully outlying)
#' @return list(beta, Gamma, gamma) mapping X -> beta * X %*% Gamma + 1 gamma'
weighted_procrustes <- function(X, Y, w) {
  k <- nrow(X); p <- ncol(X)
  w <- pmax(w, 0)
  W <- sum(w)
  if (W <= 1e-8) {
    return(list(beta = 1, Gamma = diag(p), gamma = rep(0, p)))
  }

  xc <- colSums(w * X) / W
  yc <- colSums(w * Y) / W
  Xc <- sweep(X, 2, xc, "-")
  Yc <- sweep(Y, 2, yc, "-")

  sw <- sqrt(w)
  A <- sw * Xc                      # D^{1/2} (X - xbar)
  B <- sw * Yc                      # D^{1/2} (Y - ybar)

  M  <- crossprod(A, B)             # A' B
  sv <- svd(M)
  Gamma <- sv$u %*% t(sv$v)

  # similarity transforms use proper rotations only (no reflection)
  if (p > 1 && det(Gamma) < 0) {
    Sm <- diag(c(rep(1, p - 1), -1))
    Gamma <- sv$u %*% Sm %*% t(sv$v)
  }

  num   <- sum(diag(t(Gamma) %*% M))
  den   <- sum(A^2)
  beta  <- if (den > 1e-10) num / den else 1
  gamma <- yc - beta * as.numeric(xc %*% Gamma)

  list(beta = beta, Gamma = Gamma, gamma = gamma)
}

apply_similarity <- function(X, fit) {
  sweep(fit$beta * (X %*% fit$Gamma), 2, fit$gamma, "+")
}

# ---------------------------------------------------------------------------
# 2. Weighted generalized Procrustes analysis (GPA): aligns a whole sample to
#    a common (weighted) mean shape, using an arbitrary landmark-weight list
# ---------------------------------------------------------------------------
weighted_gpa <- function(X_list, w_list, max_iter = 15, tol = 1e-8) {
  n <- length(X_list); k <- nrow(X_list[[1]]); p <- ncol(X_list[[1]])
  mu <- Reduce(`+`, X_list) / n
  mu <- mu / sqrt(sum(mu^2))

  aligned <- X_list
  for (iter in seq_len(max_iter)) {
    for (i in seq_len(n)) {
      fit <- weighted_procrustes(X_list[[i]], mu, w_list[[i]])
      aligned[[i]] <- apply_similarity(X_list[[i]], fit)
    }
    newmu <- matrix(0, k, p)
    for (l in seq_len(k)) {
      wl <- vapply(w_list, `[`, numeric(1), l)
      wsum <- sum(wl)
      if (wsum > 1e-8) {
        acc <- Reduce(`+`, Map(function(X, w) w * X[l, ], aligned, as.list(wl)))
        newmu[l, ] <- acc / wsum
      } else {
        newmu[l, ] <- mu[l, ]
      }
    }
    nrm <- sqrt(sum(newmu^2))
    if (nrm > 1e-10) newmu <- newmu / nrm
    delta <- sum((newmu - mu)^2)
    mu <- newmu
    if (delta < tol) break
  }
  list(mu = mu, aligned = aligned)
}

# ---------------------------------------------------------------------------
# 3. Doubly weighted Nadaraya-Watson shape predictor, eq (7)
# ---------------------------------------------------------------------------
gaussian_kernel <- function(u, h) dnorm(u / h)

#' @param t0 time at which to predict
#' @param times vector of observed times (length n)
#' @param X_list list of n aligned k x p configurations (the X-tilde's)
#' @param w_list list of n length-k landmark weights (delta * u)
#' @param h bandwidth
#' @param exclude index to drop (for leave-one-out evaluation), or NULL
nw_predict_shape <- function(t0, times, X_list, w_list, h, exclude = NULL) {
  k <- nrow(X_list[[1]]); p <- ncol(X_list[[1]])
  idx <- seq_along(times)
  if (!is.null(exclude)) idx <- setdiff(idx, exclude)
  Kv <- gaussian_kernel(t0 - times[idx], h)

  pred <- matrix(NA_real_, k, p)
  for (l in seq_len(k)) {
    wl  <- vapply(idx, function(i) w_list[[i]][l], numeric(1))
    wk  <- wl * Kv
    den <- sum(wk)
    if (den > 1e-10) {
      num <- Reduce(`+`, Map(function(i, g) g * X_list[[i]][l, ], idx, as.list(wk)))
      pred[l, ] <- num / den
    }
  }
  pred
}

# ---------------------------------------------------------------------------
# 4. Two-pass pilot -> Tukey biweight reweighting scheme (Section 3.2)
# ---------------------------------------------------------------------------
#' @param times observed times
#' @param X_raw list of n raw (possibly contaminated) k x p configurations
#' @param delta_list list of n length-k 0/1 missingness indicators
#' @param robust if FALSE, skips the biweight step ("aware" estimator)
#' @param h_pilot bandwidth used only for the internal pilot smooth
fit_shape_model <- function(times, X_raw, delta_list, robust = TRUE,
                             c_tukey = 4.685, h_pilot = NULL) {
  n <- length(times); k <- nrow(X_raw[[1]])
  # scale-free default: a fixed fraction of the observed time range, so
  # this works whether times live in [0,1] (simulation) or e.g. [7,150]
  # (days, as in real growth data) -- a hardcoded literal here silently
  # breaks on any other time scale (all kernel weights collapse to ~0).
  if (is.null(h_pilot)) h_pilot <- diff(range(times)) / 6

  # Step 0 / pass 1: missing-data-aware (but non-robust) alignment
  gpa0 <- weighted_gpa(X_raw, delta_list)
  Xtilde0 <- gpa0$aligned

  if (!robust) {
    return(list(Xtilde = Xtilde0, weights = delta_list, mu = gpa0$mu, sigma_hat = NA))
  }

  # Step 1: pilot smooth at every observed time (LEAVE-ONE-OUT: a
  # configuration must not be smoothed using itself, or its residual -
  # and hence sigma_hat and its own robustness weight - would be
  # optimistically shrunk toward zero; see accompanying referee report)
  pilot <- lapply(seq_len(n), function(i)
    nw_predict_shape(times[i], times, Xtilde0, delta_list, h_pilot, exclude = i))

  # Step 2: local residuals (landmark level)
  resid_list <- vector("list", n)
  all_r <- numeric(0)
  for (i in seq_len(n)) {
    e <- Xtilde0[[i]] - pilot[[i]]
    r <- sqrt(rowSums(e^2))
    r[delta_list[[i]] == 0] <- NA
    resid_list[[i]] <- r
    all_r <- c(all_r, r[!is.na(r)])
  }

  # Step 3: robust scale, calibrated to the CORRECT reference distribution.
  # NOTE (fixed per referee report): r_{i,l} = ||e_{i,l}|| is the norm of a
  # p-dimensional residual, not a scalar Gaussian deviation, so the classical
  # 1-D MAD correction 1.4826 = 1/qnorm(0.75) does NOT apply here -- using it
  # silently assumes r_{i,l} behaves like a folded-normal |Z|, which it does
  # not. If e_{i,l} ~ N_p(0, sigma^2 I), then r_{i,l}/sigma follows a chi
  # distribution with p degrees of freedom, so the consistent scale estimate
  # is median(r)/median(chi_p), not 1.4826*median(r). For p=2 this divisor is
  # sqrt(2*log(2)) =~ 1.177 (median of chi_p, not the 1-D value 0.6745), so
  # the previous formula overstated sigma_hat by a factor of about 1.75 --
  # substantially loosening the Tukey biweight cutoff c*sigma_hat away from
  # its nominal 95%-efficiency calibration (Section 3.2 / Eq. 6).
  p_dim <- ncol(X_raw[[1]])
  chi_p_median <- sqrt(stats::qchisq(0.5, df = p_dim))
  sigma_hat <- stats::median(all_r, na.rm = TRUE) / chi_p_median
  if (!is.finite(sigma_hat) || sigma_hat <= 0) sigma_hat <- 1e-6

  # Step 4: Tukey biweight down-weighting, eq (6)
  u_list <- lapply(seq_len(n), function(i) {
    r <- resid_list[[i]]
    u <- rep(0, k)
    ok <- !is.na(r)
    z <- r[ok] / (c_tukey * sigma_hat)
    uu <- (1 - z^2)^2
    uu[abs(z) > 1] <- 0
    u[ok] <- uu
    u
  })
  w_list <- Map(function(d, u) d * u, delta_list, u_list)

  # Step 5: refit alignment with the combined weight
  gpa1 <- weighted_gpa(X_raw, w_list)

  list(Xtilde = gpa1$aligned, weights = w_list, mu = gpa1$mu, sigma_hat = sigma_hat)
}

# ---------------------------------------------------------------------------
# 5. Mean-imputation for the "naive" estimator (literal original-model use)
# ---------------------------------------------------------------------------
mean_impute <- function(X_raw, delta_list) {
  n <- length(X_raw); k <- nrow(X_raw[[1]]); p <- ncol(X_raw[[1]])
  landmark_mean <- matrix(0, k, p)
  for (l in seq_len(k)) {
    obs <- t(sapply(seq_len(n), function(i) X_raw[[i]][l, ]))
    ok  <- vapply(seq_len(n), function(i) delta_list[[i]][l] == 1, logical(1))
    landmark_mean[l, ] <- if (any(ok)) colMeans(obs[ok, , drop = FALSE]) else colMeans(obs)
  }
  X_imp <- X_raw
  for (i in seq_len(n)) {
    miss <- delta_list[[i]] == 0
    if (any(miss)) X_imp[[i]][miss, ] <- landmark_mean[miss, , drop = FALSE]
  }
  X_imp
}

# ---------------------------------------------------------------------------
# 6. Hexagon growth simulation (k = 6) and contamination (Table 1)
# ---------------------------------------------------------------------------
simulate_hexagon <- function(n = 30, sigma_noise = 0.05) {
  times  <- seq(0, 1, length.out = n)
  angles <- seq(0, 5) * pi / 3
  base   <- cbind(cos(angles), sin(angles))

  X_true <- vector("list", n)
  for (i in seq_len(n)) {
    t     <- times[i]
    size  <- 1 + 0.6 * t                     # nonlinear growth in size
    rot   <- t * pi / 6                      # rotation drift
    shear <- 0.3 * t                         # genuine SHAPE change (not just size/rotation)
    Sh <- matrix(c(1, 0, shear, 1), 2, 2)
    R  <- matrix(c(cos(rot), sin(rot), -sin(rot), cos(rot)), 2, 2)
    X_true[[i]] <- size * (base %*% Sh) %*% R
  }
  X_noisy <- lapply(X_true, function(X)
    X + matrix(rnorm(length(X), 0, sigma_noise), nrow(X), ncol(X)))

  list(times = times, X_true = X_true, X_noisy = X_noisy)
}

contaminate <- function(X_list, outlier_prob = 0, missing_prob = 0, outlier_sd = 3) {
  n <- length(X_list); k <- nrow(X_list[[1]]); p <- ncol(X_list[[1]])
  X_out <- X_list
  delta_list <- replicate(n, rep(1, k), simplify = FALSE)
  for (i in seq_len(n)) {
    for (l in seq_len(k)) {
      u <- runif(1)
      if (u < missing_prob) {
        delta_list[[i]][l] <- 0
      } else if (u < missing_prob + outlier_prob) {
        X_out[[i]][l, ] <- X_out[[i]][l, ] + rnorm(p, 0, outlier_sd)
      }
    }
  }
  list(X = X_out, delta = delta_list)
}

# ---------------------------------------------------------------------------
# 7. Evaluation criteria: Procrustes RMSD, R^2, and a Peng-Deng-style
#    similarity index (Sim). All against the noise-free ground truth.
# ---------------------------------------------------------------------------
procrustes_align_full <- function(Xhat, Xtrue) {
  fit <- weighted_procrustes(Xhat, Xtrue, rep(1, nrow(Xhat)))
  apply_similarity(Xhat, fit)
}

eval_scenario <- function(pred_list, true_list, mean_true_shape) {
  n <- length(pred_list)
  ss_res <- 0; ss_tot <- 0; rmsd_v <- numeric(n); sim_v <- numeric(n)

  for (i in seq_len(n)) {
    if (anyNA(pred_list[[i]])) next
    Xa <- procrustes_align_full(pred_list[[i]], true_list[[i]])
    resid <- sum((Xa - true_list[[i]])^2)
    ss_res <- ss_res + resid

    Ma <- procrustes_align_full(mean_true_shape, true_list[[i]])
    ss_tot <- ss_tot + sum((Ma - true_list[[i]])^2)

    k <- nrow(Xa)
    rmsd_v[i] <- sqrt(resid / k)
    # Peng & Deng (2013)-style normalized similarity: 1 - relative shape
    # distance, scaled to [0,100]; a simple, monotone stand-in for the
    # original (graph-theoretic) index, since RMSD/R2 already capture the
    # essential comparison used in the paper.
    scale_fac <- sqrt(sum(sweep(true_list[[i]], 2, colMeans(true_list[[i]]))^2))
    sim_v[i] <- 100 * max(0, 1 - sqrt(resid) / scale_fac)
  }
  R2 <- 100 * (1 - ss_res / ss_tot)
  list(R2 = R2, RMSD = 100 * mean(rmsd_v), Sim = mean(sim_v))
}

# ---------------------------------------------------------------------------
# 7b. Robust summary of a (possibly heavily right-tailed) vector of replicate
# results. Under outlier contamination, R2 is bounded above by 100 but
# unbounded below, so its distribution across replications can be strongly
# skewed and a handful of catastrophic replicates can dominate the plain
# mean (see the accompanying referee report, Major point 3). We therefore
# report the median and a 10%-trimmed mean alongside mean/sd/min/max so a
# reader can see whether "mean" is being driven by a few extreme draws.
# ---------------------------------------------------------------------------
robust_summary <- function(x, trim = 0.10) {
  c(min = min(x), mean = mean(x), median = stats::median(x),
    trimmed = mean(x, trim = trim), max = max(x), sd = sd(x))
}

# ---------------------------------------------------------------------------
# 7c. Paired significance test between two methods evaluated with COMMON
# RANDOM NUMBERS (i.e. the same simulated/contaminated draw in each
# replicate; see the "common random numbers" restructuring in
# run_simulation()/run_real_data_contamination_study() below). Because the
# two vectors are matched replicate-for-replicate, a paired Wilcoxon
# signed-rank test is valid and much more powerful than an unpaired
# comparison of means (see referee report, Major point 5). Returns NA if
# fewer than 2 finite paired differences are available (e.g. identical
# vectors, or all-NA results).
# ---------------------------------------------------------------------------
paired_method_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2 || all(x[ok] == y[ok])) return(NA_real_)
  # If both series are individually constant across "replications" (no
  # random contamination was actually drawn -- e.g. a Clean scenario that
  # reuses the same fixed, uncontaminated input every replicate), the
  # n_rep repeats are not independent random draws. A paired test on them
  # is not a valid statistical inference even though x and y differ from
  # each other by a fixed nonzero amount; report NA rather than a
  # degenerate p-value (see referee report).
  if (stats::sd(x[ok]) == 0 && stats::sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::wilcox.test(x[ok], y[ok], paired = TRUE)$p.value)
}

# ---------------------------------------------------------------------------
# 8. One leave-one-out replication for a given scenario / method
# ---------------------------------------------------------------------------
run_method_once <- function(times, X_raw, delta_list, X_true, method,
                             h_grid = seq(0.02, 0.45, by = 0.03)) {
  n <- length(times)

  if (method == "naive") {
    X_used <- mean_impute(X_raw, delta_list)
    delta_used <- replicate(n, rep(1, nrow(X_raw[[1]])), simplify = FALSE)
    fitted <- fit_shape_model(times, X_used, delta_used, robust = FALSE)
  } else if (method == "aware") {
    fitted <- fit_shape_model(times, X_raw, delta_list, robust = FALSE)
  } else if (method == "robust") {
    fitted <- fit_shape_model(times, X_raw, delta_list, robust = TRUE)
  } else stop("unknown method")

  Xtilde <- fitted$Xtilde
  W      <- fitted$weights

  # Bandwidth by leave-one-out grid search using the RLSCV criterion of
  # eq (8): a MEDIAN (not mean) over configurations of the squared
  # distance to the LOO prediction, so that a handful of still-heavily
  # contaminated points cannot dominate the bandwidth choice itself.
  # Each configuration's distance is computed only over its reliably
  # weighted landmarks, so the criterion is consistent with (rather than
  # blind to) the robustness/missingness weights used for the fit.
  best_h <- h_grid[1]; best_err <- Inf
  for (h in h_grid) {
    errs <- numeric(0)
    for (i in seq_len(n)) {
      p_i <- nw_predict_shape(times[i], times, Xtilde, W, h, exclude = i)
      wl  <- W[[i]]
      if (!anyNA(p_i) && sum(wl) > 1e-8) {
        d2 <- sum(wl * rowSums((p_i - Xtilde[[i]])^2)) / sum(wl)
        errs <- c(errs, d2)
      }
    }
    if (length(errs) > 0) {
      err <- stats::median(errs)
      if (err < best_err) { best_err <- err; best_h <- h }
    }
  }

  # leave-one-out predictions evaluated against noise-free ground truth
  preds <- lapply(seq_len(n), function(i)
    nw_predict_shape(times[i], times, Xtilde, W, best_h, exclude = i))

  mean_true <- Reduce(`+`, X_true) / n
  eval_scenario(preds, X_true, mean_true)
}

# ---------------------------------------------------------------------------
# 9. Simulation driver: reproduces the structure of Table 1 / Table 2
# ---------------------------------------------------------------------------
run_simulation <- function(n_rep = 40, n = 30, sigma_noise = 0.05, seed = 1) {
  set.seed(seed)
  scenarios <- list(
    Clean         = c(outlier_prob = 0,    missing_prob = 0),
    OutliersOnly  = c(outlier_prob = 0.08, missing_prob = 0),
    MissingOnly   = c(outlier_prob = 0,    missing_prob = 0.15),
    Both          = c(outlier_prob = 0.08, missing_prob = 0.15)
  )
  methods_all <- c("naive", "aware", "robust")

  results <- list()
  for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]
    # naive/aware are identical when there is no missingness; skip the
    # duplicate the way Table 2 of the accompanying paper does
    methods <- if (sc_name %in% c("Clean", "OutliersOnly")) c("naive", "robust") else methods_all

    # COMMON RANDOM NUMBERS across methods within a replicate (same hexagon
    # draw + same contamination draw for naive/aware/robust) -- see the
    # identical restructuring, and its rationale, in
    # run_real_data_contamination_study() above and referee report Major
    # point 5. This also makes the paired significance test below valid.
    R2m <- RMSDm <- Simm <- matrix(NA_real_, n_rep, length(methods),
                                    dimnames = list(NULL, methods))
    for (rep in seq_len(n_rep)) {
      sim  <- simulate_hexagon(n = n, sigma_noise = sigma_noise)
      cont <- contaminate(sim$X_noisy,
                           outlier_prob = sc["outlier_prob"],
                           missing_prob = sc["missing_prob"])
      for (method in methods) {
        out <- run_method_once(sim$times, cont$X, cont$delta, sim$X_true, method)
        R2m[rep, method] <- out$R2; RMSDm[rep, method] <- out$RMSD; Simm[rep, method] <- out$Sim
      }
    }

    for (method in methods) {
      rs_R2   <- robust_summary(R2m[, method])
      rs_RMSD <- robust_summary(RMSDm[, method])
      rs_Sim  <- robust_summary(Simm[, method])
      p_vs_naive <- if (method != "naive" && "naive" %in% methods) paired_method_test(R2m[, method], R2m[, "naive"]) else NA_real_
      p_vs_aware <- if (method != "aware" && "aware" %in% methods) paired_method_test(R2m[, method], R2m[, "aware"]) else NA_real_
      key <- paste(sc_name, method, sep = " / ")
      results[[key]] <- data.frame(
        Scenario = sc_name, Method = method,
        Sim_min = rs_Sim["min"],  Sim_mean = rs_Sim["mean"], Sim_median = rs_Sim["median"], Sim_max = rs_Sim["max"], Sim_sd = rs_Sim["sd"],
        RMSD_min = rs_RMSD["min"], RMSD_mean = rs_RMSD["mean"], RMSD_median = rs_RMSD["median"], RMSD_max = rs_RMSD["max"], RMSD_sd = rs_RMSD["sd"],
        R2_min = rs_R2["min"], R2_mean = rs_R2["mean"], R2_median = rs_R2["median"], R2_trimmed10 = rs_R2["trimmed"], R2_max = rs_R2["max"], R2_sd = rs_R2["sd"],
        p_R2_vs_naive = p_vs_naive, p_R2_vs_aware = p_vs_aware,
        row.names = NULL
      )
      cat(sprintf("%-25s  R2 mean=%7.2f median=%7.2f  RMSD mean=%6.2f  p(vs naive)=%s\n",
                   key, rs_R2["mean"], rs_R2["median"], rs_RMSD["mean"],
                   if (is.na(p_vs_naive)) "--" else sprintf("%.3f", p_vs_naive)))
    }
  }
  do.call(rbind, results)
}

# ---------------------------------------------------------------------------
# 10. Real-data support
# ---------------------------------------------------------------------------
# 10a. Flexible loader: turn a k x p x n array (the usual "shapes"-package
# convention, e.g. shapes::rats$x for one rat, dim 8 x 2 x 8) into the
# list-of-matrices format used throughout this script. Also accepts a
# k x p x n x N array (N repeated units, e.g. all 18 rats) and returns a
# list of such n-length lists, one per unit.
array_to_configs <- function(arr) {
  d <- dim(arr)
  if (length(d) == 3) {
    k <- d[1]; p <- d[2]; n <- d[3]
    lapply(seq_len(n), function(i) arr[, , i])
  } else if (length(d) == 4) {
    k <- d[1]; p <- d[2]; n <- d[3]; N <- d[4]
    lapply(seq_len(N), function(u)
      lapply(seq_len(n), function(i) arr[, , i, u]))
  } else stop("expected a k x p x n or k x p x n x N array")
}

# 10b. Direct application to real (assumed clean) data: fits naive / aware /
# robust and reports in-sample and leave-one-out criteria, in the spirit of
# Sections 3.2-3.3 of the original Moghimbeygi & Golalizadeh (2024) paper.
# Use this when the real configurations have NO known missingness (delta_list
# defaults to all-observed) - e.g. to reproduce their rat-skull / DNA analysis
# with the (here, un-robustified) "aware" method reducing to their original
# estimator, and to see what the robust method changes on real, un-contaminated
# data (it should change very little, per the "no cost when clean" property).
fit_real_data <- function(times, X_list, delta_list = NULL, method = "robust",
                           h_grid = NULL) {
  n <- length(times); k <- nrow(X_list[[1]])
  if (is.null(delta_list)) delta_list <- replicate(n, rep(1, k), simplify = FALSE)
  if (is.null(h_grid)) {
    span <- diff(range(times))
    h_grid <- seq(span / 40, span / 2, length.out = 15)
  }
  out <- run_method_once(times, X_list, delta_list, X_list, method, h_grid = h_grid)
  out
}

# 10c. Contamination study on REAL data: treats the real, observed
# configurations as ground truth and superimposes the SAME synthetic
# missingness / outlier scheme as Table 1 on top of them, then compares
# naive / aware / robust exactly as in run_simulation() above. This is the
# validation the original short paper deferred to "future work": it
# directly checks whether the method's simulated gains carry over to a
# real morphometric / dynamical dataset's own landmark configuration and
# growth trajectory, rather than only the toy hexagon.
#
# 10c-i. Explicit, documented outlier-magnitude formula for real data.
# There is no known noise level for real digitized landmarks (unlike the
# hexagon simulation, which has a known sigma_noise), so we anchor the
# outlier standard deviation to a multiple of the data's own typical
# landmark-to-centroid spread:
#     typical_scale = mean over configurations i of
#                      sqrt( sum_l || X_i(l,.) - mean_l X_i(l,.) ||^2 / k )
#     outlier_sd     = outlier_sd_factor * typical_scale
# With outlier_sd_factor = 3 this places an outlier landmark, on average,
# about 3 "configuration radii" away from where it should be -- the same
# proportion (outlier sd / object scale) used for the hexagon, where
# outlier sd = 3 against a hexagon of circumradius ~= 1. This function is
# the single source of truth for that formula (an earlier version of this
# script computed an intermediate `outlier_sd` value that was never
# actually used -- see the accompanying referee report, Major point 2 --
# that dead code has been removed).
contamination_scale <- function(X_list, outlier_sd_factor = 3) {
  typical_scale <- mean(vapply(X_list, function(X)
    sqrt(sum(sweep(X, 2, colMeans(X))^2) / nrow(X)), numeric(1)))
  list(typical_scale = typical_scale, outlier_sd = outlier_sd_factor * typical_scale)
}

# @param times   observed time vector (length n)
# @param X_list  list of n k x p REAL configurations (ground truth)
run_real_data_contamination_study <- function(times, X_list, n_rep = 40,
                                               outlier_sd_factor = 3, seed = 1) {
  set.seed(seed)
  n <- length(times); k <- nrow(X_list[[1]])

  cs <- contamination_scale(X_list, outlier_sd_factor)
  typical_scale <- cs$typical_scale
  outlier_sd    <- cs$outlier_sd

  scenarios <- list(
    Clean         = c(outlier_prob = 0,    missing_prob = 0),
    OutliersOnly  = c(outlier_prob = 0.08, missing_prob = 0),
    MissingOnly   = c(outlier_prob = 0,    missing_prob = 0.15),
    Both          = c(outlier_prob = 0.08, missing_prob = 0.15)
  )
  methods_all <- c("naive", "aware", "robust")
  span <- diff(range(times))
  h_grid <- seq(span / 40, span / 2, length.out = 12)

  results <- list()
  for (sc_name in names(scenarios)) {
    sc <- scenarios[[sc_name]]
    methods <- if (sc_name %in% c("Clean", "OutliersOnly")) c("naive", "robust") else methods_all

    # COMMON RANDOM NUMBERS: draw ONE contaminated dataset per replicate and
    # evaluate every method of this scenario on that SAME draw, rather than
    # re-randomizing independently per method. This (i) reduces Monte Carlo
    # variance in the naive-vs-robust comparison and (ii) is what makes the
    # paired significance test below valid (see referee report, Major
    # point 5; the previous version of this function re-contaminated the
    # data independently inside the per-method loop, so replicate #5 of
    # "naive" and replicate #5 of "robust" were not comparable pairs).
    R2m <- RMSDm <- Simm <- matrix(NA_real_, n_rep, length(methods),
                                    dimnames = list(NULL, methods))
    for (rep in seq_len(n_rep)) {
      cont <- contaminate(X_list, outlier_prob = sc["outlier_prob"],
                           missing_prob = sc["missing_prob"],
                           outlier_sd = outlier_sd)
      for (method in methods) {
        out <- run_method_once(times, cont$X, cont$delta, X_list, method, h_grid = h_grid)
        R2m[rep, method] <- out$R2; RMSDm[rep, method] <- out$RMSD; Simm[rep, method] <- out$Sim
      }
    }

    for (method in methods) {
      rs_R2   <- robust_summary(R2m[, method])
      rs_RMSD <- robust_summary(RMSDm[, method])
      rs_Sim  <- robust_summary(Simm[, method])
      p_vs_naive  <- if (method != "naive"  && "naive"  %in% methods) paired_method_test(R2m[, method], R2m[, "naive"])  else NA_real_
      p_vs_aware  <- if (method != "aware"  && "aware"  %in% methods) paired_method_test(R2m[, method], R2m[, "aware"])  else NA_real_
      key <- paste(sc_name, method, sep = " / ")
      results[[key]] <- data.frame(
        Scenario = sc_name, Method = method,
        Sim_mean = rs_Sim["mean"], Sim_median = rs_Sim["median"],
        RMSD_mean = rs_RMSD["mean"], RMSD_median = rs_RMSD["median"],
        R2_min = rs_R2["min"], R2_mean = rs_R2["mean"], R2_median = rs_R2["median"],
        R2_trimmed10 = rs_R2["trimmed"], R2_max = rs_R2["max"], R2_sd = rs_R2["sd"],
        p_R2_vs_naive = p_vs_naive, p_R2_vs_aware = p_vs_aware,
        row.names = NULL
      )
      cat(sprintf("%-25s  R2 mean=%7.2f median=%7.2f  RMSD mean=%6.2f  p(vs naive)=%s\n",
                   key, rs_R2["mean"], rs_R2["median"], rs_RMSD["mean"],
                   if (is.na(p_vs_naive)) "--" else sprintf("%.3f", p_vs_naive)))
    }
  }
  do.call(rbind, results)
}

# ---------------------------------------------------------------------------
# Usage on your own real data
# ---------------------------------------------------------------------------
# Option A - the "shapes" R package (Dryden, 2019), as used in the original
# paper (install.packages("shapes"); this environment has no CRAN access,
# so this must be run on your own machine):
#
#   library(shapes)
#   data(rats)                       # rats$x: 8 landmarks x 2 dims x 8 ages x 18 rats
#   arr    <- rats$x[, , , 2]        # rat number 2 -> 8 x 2 x 8 array
#   times  <- c(7, 14, 21, 30, 40, 60, 90, 150)
#   X_rat2 <- array_to_configs(arr)  # list of 8 (8 x 2) configurations
#   fit_real_data(times, X_rat2, method = "robust")
#   study  <- run_real_data_contamination_study(times, X_rat2, n_rep = 40)
#
# Option B - your own data, already as an RDS/csv you can turn into a list
# of k x p matrices in R (one matrix per observed time), plus a numeric
# `times` vector of the same length:
#
#   X_list <- readRDS("my_configurations.rds")   # list of k x p matrices
#   times  <- readRDS("my_times.rds")
#   fit_real_data(times, X_list, method = "robust")
#   study  <- run_real_data_contamination_study(times, X_list, n_rep = 40)
#
# Both fit_real_data() and run_real_data_contamination_study() work with ANY
# k (landmark count) and p (dimension: 2 for planar, 3 for DNA/skull data),
# since the weighted-SVD solution (5) does not depend on the ambient dimension.

# ---------------------------------------------------------------------------
# Run a (small, fast) demonstration if executed as a script directly.
# For a full Table-2-scale run, call run_simulation(n_rep = 40).
# ---------------------------------------------------------------------------
if (sys.nframe() == 0) {
  tab <- run_simulation(n_rep = 10)   # small n_rep for a quick demo run
  print(tab, row.names = FALSE)
}

