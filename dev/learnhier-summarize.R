# Lane `learnhier` (item 2.2): turn a directory of replicate records
# into the recovery table.
#
#   Rscript dev/learnhier-summarize.R <dir> <design> [<design> ...]
#
# TWO TRUTHS PER COMPONENT, and the difference between them is the
# finding this row exists to record.
#
#   population  the value the design was written with. For every
#               component of the bandit block and for the three rlddm
#               components whose fitted link IS the drawn link, this is
#               a truth in the ordinary sense and the coverage column
#               beside it is a coverage.
#   realized    the same component computed from THIS replicate's own
#               drawn deviations and its own observed floors, with no
#               fit. For `ndt` under ndt_group() there is no population
#               value to have: the link's ceiling is each learner's
#               fastest response, which is a property of the data. The
#               realized column is what the model is actually estimating
#               there.
#
# Coverage against a realized target is NOT a coverage and is not read
# as one: the target carries the sampling noise of 100 draws and, in the
# limit where the deviations were observed rather than estimated, the
# estimate would equal it exactly and the rate would go to one. It is
# reported because a rate far BELOW nominal there still says the fit is
# not reading the block it fitted.

source("dev/learnhier-env.R")
args <- commandArgs(trailingOnly = TRUE)
DIR <- args[[1L]]
DESIGNS <- args[-1L]

# The population truths, by design and by component name. A name absent
# here has no population value and is scored against the realized column
# alone.
lh_pop_truth <- function(design) {
  if (design %in% c("bandit", "bandit0")) {
    cr <- if (design == "bandit") 0.5 else 0
    c("b_alpha.(Intercept)" = stats::qlogis(0.35),
      "b_tau.(Intercept)" = log(3),
      sd_alpha = 0.5, sd_tau = 0.3, "cor_alpha~tau" = cr)
  } else {
    c("b_alpha.(Intercept)" = stats::qlogis(0.35),
      "b_drift.(Intercept)" = 2.5,
      "b_bs.(Intercept)" = log(1.5),
      sd_alpha = 0.5, sd_drift = 1.0, sd_bs = 0.2,
      "cor_alpha~drift" = 0.4, "cor_alpha~bs" = 0.0,
      "cor_drift~bs" = -0.3)
  }
}

# The realized target for each component, read off one record. The
# fixed-effect intercepts have no realized counterpart except `ndt`'s,
# whose oracle value IS the realized one.
lh_realized <- function(r) {
  out <- r$truth_fitted
  if (identical(r$design, "rlddm") && length(r$oracle_betad) >= 3L) {
    out[["b_ndt.(Intercept)"]] <- r$oracle_betad[[3L]]
  }
  out
}

lh_read <- function(dir, design) {
  fs <- list.files(dir, pattern = paste0("^", design, "-[0-9]+[.]rds$"),
                   full.names = TRUE)
  lapply(fs, readRDS)
}

# One component's row. `est`, `lwr`, `upr` come from the fit; `pop` and
# `real` are the two targets; every spread is stated as something this
# run measured.
lh_row <- function(name, est, lwr, upr, pop, real) {
  n <- sum(is.finite(est))
  sdv <- stats::sd(est, na.rm = TRUE)
  # a component on the boundary of the parameter space gets no interval
  # at all, and a table that silently dropped those replicates would
  # report a coverage over the ones that behaved
  n_no_ci <- sum(!(is.finite(lwr) & is.finite(upr)))
  cov_of <- function(tg) {
    ok <- is.finite(lwr) & is.finite(upr) & is.finite(tg)
    if (!any(ok)) return(c(NA_real_, NA_real_, 0))
    hit <- lwr[ok] <= tg[ok] & upr[ok] >= tg[ok]
    p <- mean(hit)
    c(p, sqrt(p * (1 - p) / sum(ok)), sum(ok))
  }
  cp <- cov_of(pop)
  cr <- cov_of(real)
  data.frame(component = name, n = n, no_ci = n_no_ci,
             pop_truth = if (all(is.na(pop))) NA_real_ else pop[[1L]],
             mean_est = mean(est, na.rm = TRUE),
             min_est = min(est, na.rm = TRUE),
             max_est = max(est, na.rm = TRUE),
             bias_pop = mean(est - pop, na.rm = TRUE),
             mcse = sdv / sqrt(n), sd_est = sdv,
             cov_pop = cp[[1L]], cov_pop_se = cp[[2L]],
             mean_real = mean(real, na.rm = TRUE),
             bias_real = mean(est - real, na.rm = TRUE),
             cov_real = cr[[1L]], cov_real_se = cr[[2L]],
             stringsAsFactors = FALSE)
}

lh_summary <- function(recs, design) {
  ok <- vapply(recs, function(r) isTRUE(r$ok), TRUE)
  cat("\n================ ", design, " ================\n", sep = "")
  cat("records ", length(recs), ", fitted ", sum(ok), "\n", sep = "")
  if (!any(ok)) {
    for (r in recs[!ok]) cat("  FAILED seed ", r$seed, ": ",
                             substr(r$why, 1L, 160L), "\n", sep = "")
    return(invisible(NULL))
  }
  for (r in recs[!ok]) cat("  FAILED seed ", r$seed, ": ",
                           substr(r$why, 1L, 160L), "\n", sep = "")
  R <- recs[ok]
  g <- function(f) vapply(R, f, numeric(1))
  cat(R[[1L]]$ns, " learners x ", R[[1L]]$nt, " trials, ",
      R[[1L]]$rows, " rows, ", R[[1L]]$n_par,
      " outer parameters\n", sep = "")
  # EVERY REPLICATE COUNTED HERE IS A FILE ON DISK. A launcher's own
  # tally is not a replicate count: this lane lost a 44-minute rlddm fit
  # that its launcher had counted as started and that wrote nothing, so
  # the audit below states the distinct seeds actually present and the
  # gaps against the seed grid the launcher was given.
  sd_all <- g(function(r) r$seed)
  grid <- seq(min(sd_all), max(sd_all), by = 10L)
  cat("seeds: ", length(unique(sd_all)), " distinct of ", length(sd_all),
      " files, ", min(sd_all), " to ", max(sd_all), " step 10",
      if (length(setdiff(grid, sd_all))) {
        paste0("; MISSING ", length(setdiff(grid, sd_all)), ": ",
               paste(utils::head(setdiff(grid, sd_all), 8L),
                     collapse = ","))
      } else "; no gaps", "\n", sep = "")
  secs <- g(function(r) r$fit_s)
  cat("fit seconds: min ", round(min(secs), 1), ", median ",
      round(stats::median(secs), 1), ", max ", round(max(secs), 1),
      ", total ", round(sum(secs) / 60, 1), " min\n", sep = "")
  cat("convergence code 0: ", sum(g(function(r) r$conv) == 0L), " of ",
      length(R), "; positive definite Hessian: ",
      sum(vapply(R, function(r) isTRUE(r$pdHess), TRUE)), "; NaN se: ",
      sum(g(function(r) r$n_bad_se)), "\n", sep = "")
  mg <- g(function(r) r$max_grad)
  cat("max |gradient|: median ", format(stats::median(mg), digits = 3),
      ", worst ", format(max(mg), digits = 3), "\n", sep = "")
  dll <- g(function(r) r$logLik - r$oracle)
  dll <- dll[is.finite(dll)]
  if (!length(dll)) {
    cat("no oracle stated for this design\n")
  } else {
    cat("logLik minus the oracle at this replicate's own truths: min ",
        round(min(dll), 4), ", median ", round(stats::median(dll), 3),
        ", max ", round(max(dll), 2), "; below zero on ",
        sum(dll < 0), " of ", length(dll), "\n", sep = "")
  }

  nms <- unique(c(sub("[.](est|lwr|upr)$", "", names(R[[1L]]$fixed)),
                  sub("[.](est|lwr|upr)$", "", names(R[[1L]]$vc))))
  pop <- lh_pop_truth(design)
  rows <- list()
  att <- list()
  for (nm in nms) {
    src <- if (nm %in% sub("[.](est|lwr|upr)$", "",
                           names(R[[1L]]$fixed))) "fixed" else "vc"
    pull <- function(sfx) {
      key <- paste0(nm, ".", sfx)
      vapply(R, function(r) {
        v <- r[[src]]
        if (key %in% names(v)) v[[key]] else NA_real_
      }, numeric(1))
    }
    est <- pull("est")
    if (all(!is.finite(est))) next
    if (all(est == 0)) next
    real <- vapply(R, function(r) {
      v <- lh_realized(r)
      if (nm %in% names(v)) v[[nm]] else NA_real_
    }, numeric(1))
    pv <- if (nm %in% names(pop)) rep(pop[[nm]], length(R)) else {
      rep(NA_real_, length(R))
    }
    rows[[nm]] <- lh_row(nm, est, pull("lwr"), pull("upr"), pv, real)
    att[[nm]] <- list(est = est, lo = pull("lwr"), hi = pull("upr"),
                      pop = pv)
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  num <- vapply(tab, is.numeric, TRUE)
  tab[num] <- lapply(tab[num], function(x) round(x, 4))
  cat("\n")
  print(tab[c("component", "n", "no_ci", "pop_truth", "mean_est",
              "bias_pop", "mcse", "sd_est", "cov_pop", "cov_pop_se")])
  cat("\nthe same components against the REALIZED target",
      " (see the header of this script for why that is not a coverage)\n",
      sep = "")
  print(tab[c("component", "mean_real", "bias_real", "cov_real",
              "cov_real_se", "min_est", "max_est")])
  # EXCLUSION ON THE FLAGS, because every coverage above is built from
  # Wald intervals and those are built from standard errors a
  # non-positive-definite Hessian does not license. The convergence code
  # and the Hessian check are NOT the same set, so both are reported.
  flags <- list(
    "convergence code 0" = vapply(R, function(r) r$conv == 0L, TRUE),
    "positive definite Hessian" = vapply(R, function(r) isTRUE(r$pdHess),
                                         TRUE))
  for (nm in names(flags)) {
    keep <- flags[[nm]]
    if (all(keep) || sum(keep) < 3L) {
      cat("
exclusion on ", nm, ": nothing to exclude
", sep = "")
      next
    }
    sub <- do.call(rbind, lapply(names(rows), function(k) {
      r0 <- rows[[k]]
      est <- att[[k]]$est[keep]
      lo <- att[[k]]$lo[keep]
      hi <- att[[k]]$hi[keep]
      pv <- att[[k]]$pop[keep]
      ok <- is.finite(lo) & is.finite(hi) & is.finite(pv)
      cv <- if (any(ok)) mean(lo[ok] <= pv[ok] & hi[ok] >= pv[ok]) else NA
      data.frame(component = k, cov_all = r0$cov_pop, cov_kept = cv,
                 stringsAsFactors = FALSE)
    }))
    sub$delta <- sub$cov_kept - sub$cov_all
    cat("
exclusion on ", nm, ": dropping ", sum(!keep), " of ",
        length(keep), " replicates
", sep = "")
    cat("  largest change in any coverage against a population truth: ",
        format(max(abs(sub$delta), na.rm = TRUE), digits = 3),
        "; components whose coverage moves by more than the ",
        "binomial standard error at this count (",
        format(sqrt(0.95 * 0.05 / sum(keep)), digits = 2), "): ",
        sum(abs(sub$delta) > sqrt(0.95 * 0.05 / sum(keep)),
            na.rm = TRUE), "
", sep = "")
    print(sub[!is.na(sub$cov_kept), ], row.names = FALSE)
  }

  out <- Sys.getenv("LEARNHIER_TSV", "")
  if (nzchar(out)) {
    utils::write.table(cbind(design = design, tab),
                       file.path(out, paste0("learnhier-", design, ".tsv")),
                       sep = "\t", row.names = FALSE, quote = FALSE)
  }

  if (!is.null(R[[1L]]$ndt)) {
    nd <- do.call(rbind, lapply(R, function(r) r$ndt))
    cat("\nper-learner non-decision time, at ", R[[1L]]$nt,
        " trials (a floor is an OBSERVED minimum and moves with this",
        " count)\n", sep = "")
    q <- function(x) c(mean = mean(x), sd = stats::sd(x),
                       min = min(x), max = max(x))
    st <- t(apply(nd, 2L, q))
    print(round(st, 4))
  }
  invisible(tab)
}

for (dz in DESIGNS) {
  recs <- lh_read(DIR, dz)
  if (!length(recs)) {
    cat("\nno records for ", dz, " in ", DIR, "\n", sep = "")
    next
  }
  lh_summary(recs, dz)
}
