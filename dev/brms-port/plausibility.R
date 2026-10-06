# Estimate plausibility: frmtmb's maximum-likelihood fit of every model
# that runs, against brms's posterior on the same data.
#
#   PORT_OUT=<dir> Rscript plausibility.R
#
# Reads <PORT_OUT>/results-spell (run-vignette.R, frmtmb) and
# <PORT_OUT>/brms (brms-fit.R). Both runners seed every expression from
# its id (shim.R), so a vignette that simulates its data gives both the
# same rows.
#
# The yardstick is the brms posterior SD of each parameter: for a
# coefficient, z = (frmtmb estimate - brms posterior mean) / brms
# posterior SD. ML and a posterior mean differ by prior and by
# skewness, so |z| well under 1 is agreement; a boundary variance
# component, where a prior holds the posterior off zero, is expected
# to sit far out and is reported in its own column.
HERE <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
source(file.path(HERE, "env.R"))
VIGS <- c("brms_overview", "brms_multilevel", "brms_distreg", "brms_nonlinear",
          "brms_phylogenetics", "brms_monotonic", "brms_multivariate",
          "brms_missings", "brms_customfamilies")
rd <- function(d) {
  out <- list()
  for (v in VIGS) {
    f <- file.path(PORT_OUT, d, paste0(v, ".rds"))
    if (file.exists(f)) out <- c(out, readRDS(f))
  }
  out
}
FR <- rd("results-spell")
BR <- rd("brms")
# A fit refitted at the vignette's own budget (brms-fit.R with iter 0)
# replaces the cheap-budget fit of the same id.
BV <- rd("brms-vig")
for (id in names(BV)) {
  if (identical(BV[[id]]$status, "OK") && !is.null(BV[[id]]$fit)) {
    BR[[id]] <- BV[[id]]
    BR[[id]]$budget <- "vignette"
  }
}
hdr <- readLines(file.path(PORT_OUT, "brms", paste0(VIGS[1], ".log")), n = 1)
cat("## provenance\n\nfrmtmb:", sub("^# ", "", readLines(file.path(
  PORT_OUT, "results-spell", paste0(VIGS[1], ".log")), n = 1)), "\n")
cat("brms:  ", sub("^# ", "", hdr), "\n")
cat("refitted at the vignette's own budget:",
    paste(names(BR)[vapply(BR, function(r) identical(r$budget, "vignette"),
                           NA)], collapse = ", "), "\n\n")

cmp <- function(a, b, bsd) {
  k <- intersect(names(a), names(b))
  if (!length(k)) return(NULL)
  z <- (a[k] - b[k]) / bsd[k]
  data.frame(name = k, frmtmb = unname(a[k]), brms = unname(b[k]),
             brms_sd = unname(bsd[k]), z = unname(z),
             stringsAsFactors = FALSE)
}
f3 <- function(x) formatC(x, digits = 3, format = "fg", flag = "#")

rows <- list()
detail <- list()
model_ids <- names(FR)[vapply(FR, function(r) identical(r$kind, "model"), NA)]
for (id in model_ids) {
  fr <- FR[[id]]
  br <- BR[[id]]
  b_status <- if (is.null(br)) "NOT RUN" else br$status
  f_status <- fr$status
  fx <- if (!is.null(fr$fit$fixef) && !is.null(br$fit$fixef))
    cmp(fr$fit$fixef, br$fit$fixef, br$fit$fixef_sd) else NULL
  sd <- if (!is.null(fr$fit$sds) && !is.null(br$fit$sds))
    cmp(fr$fit$sds, br$fit$sds, br$fit$sds_sd) else NULL
  only_f <- setdiff(names(fr$fit$fixef), names(br$fit$fixef))
  only_b <- setdiff(names(br$fit$fixef), names(fr$fit$fixef))
  worst <- if (!is.null(fx)) fx[which.max(abs(fx$z)), ] else NULL
  rows[[id]] <- data.frame(
    id = id, frmtmb = f_status, brms = b_status,
    n_b = if (is.null(fx)) 0L else nrow(fx),
    med_z = if (is.null(fx)) NA else stats::median(abs(fx$z)),
    max_z = if (is.null(fx)) NA else max(abs(fx$z)),
    worst = if (is.null(worst)) "" else worst$name,
    n_sd = if (is.null(sd)) 0L else nrow(sd),
    max_z_sd = if (is.null(sd)) NA else max(abs(sd$z)),
    unmatched = length(only_f) + length(only_b),
    rhat = if (is.null(br$fit)) NA else br$fit$rhat_max,
    div = if (is.null(br$fit)) NA else br$fit$divergent,
    ndraws = if (is.null(br$fit)) NA else br$fit$ndraws,
    # frmtmb's own optimizer code: neither side of a comparison is a
    # reference value when the ML fit did not converge
    fconv = if (is.null(fr$fit$conv)) NA else fr$fit$conv,
    fwarn = {
      w <- grep("did not report convergence", fr$warnings, value = TRUE)
      if (length(w)) substr(w[1], 1, 90) else ""
    },
    stringsAsFactors = FALSE)
  detail[[id]] <- list(fx = fx, sd = sd, only_f = only_f, only_b = only_b,
                       f_msg = fr$msg, b_msg = br$msg)
}
R <- do.call(rbind, rows)

cat("## per model (generated)\n\n")
cat("| model | frmtmb | frmtmb conv | brms | coefs | median abs z |",
    "max abs z (coef) | sds | max abs z, sd | unmatched names |",
    "brms max Rhat | divergent | draws |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(R))) {
  r <- R[i, ]
  cat(sprintf(paste0("| `%s` | %s | %s | %s | %d | %s | %s | %d | %s |",
                     " %d | %s | %s | %s |\n"),
              sub("^brms_", "", r$id), r$frmtmb,
              if (is.na(r$fconv)) "" else r$fconv, r$brms, r$n_b,
              if (is.na(r$med_z)) "" else f3(r$med_z),
              if (is.na(r$max_z)) "" else
                paste0(f3(r$max_z), " (`", r$worst, "`)"),
              r$n_sd, if (is.na(r$max_z_sd)) "" else f3(r$max_z_sd),
              r$unmatched, if (is.na(r$rhat)) "" else f3(r$rhat),
              if (is.na(r$div)) "" else r$div,
              if (is.na(r$ndraws)) "" else r$ndraws))
}
both <- R[R$frmtmb == "OK" & R$brms == "OK" & R$n_b > 0, ]
cat(sprintf(paste0("\nmodels fitted by both with matched coefficients: %d ",
                   "of %d model calls\n"), nrow(both), nrow(R)))
allz <- unlist(lapply(detail[both$id], function(d) abs(d$fx$z)))
cat(sprintf("matched coefficients: %d; abs z: median %s, 90th pct %s, max %s\n",
            length(allz), f3(stats::median(allz)),
            f3(stats::quantile(allz, 0.9)), f3(max(allz))))
cat(sprintf("coefficients with abs z > 1: %d; > 2: %d\n",
            sum(allz > 1), sum(allz > 2)))
cat(sprintf("brms fits with max Rhat > 1.05: %d (%s); with divergences: %d\n",
            sum(both$rhat > 1.05, na.rm = TRUE),
            paste(both$id[which(both$rhat > 1.05)], collapse = ", "),
            sum(both$div > 0, na.rm = TRUE)))
# The same summary over the brms fits that converged, since a posterior
# mean from unmixed chains is not a reference value
ok <- both$id[!is.na(both$rhat) & both$rhat <= 1.05]
cz <- unlist(lapply(detail[ok], function(d) abs(d$fx$z)))
cat(sprintf(paste0("converged brms fits only (max Rhat <= 1.05): %d models, ",
                   "%d coefficients; abs z median %s, 90th pct %s, max %s; ",
                   "> 1: %d\n"),
            length(ok), length(cz), f3(stats::median(cz)),
            f3(stats::quantile(cz, 0.9)), f3(max(cz)), sum(cz > 1)))
nf <- both[!is.na(both$fconv) & both$fconv != 0, ]
cat(sprintf("frmtmb fits that did not converge (opt$convergence != 0): %d\n",
            nrow(nf)))
for (i in seq_len(nrow(nf))) {
  cat(sprintf("  %s: code %s, %s\n", nf$id[i], nf$fconv[i], nf$fwarn[i]))
}
# Both sides converged: brms max Rhat <= 1.05 and frmtmb code 0.
ok2 <- both$id[!is.na(both$rhat) & both$rhat <= 1.05 &
                 !is.na(both$fconv) & both$fconv == 0]
cz2 <- unlist(lapply(detail[ok2], function(d) abs(d$fx$z)))
cat(sprintf(paste0("both converged: %d models, %d coefficients; abs z median ",
                   "%s, 90th pct %s, max %s; > 1: %d\n"),
            length(ok2), length(cz2), f3(stats::median(cz2)),
            f3(stats::quantile(cz2, 0.9)), f3(max(cz2)), sum(cz2 > 1)))
sds_of <- function(ids) {
  do.call(rbind, lapply(ids, function(id) {
    d <- detail[[id]]$sd
    if (is.null(d)) NULL else cbind(id = id, d)
  }))
}
for (set in list(list("brms converged", ok), list("both converged", ok2))) {
  S <- sds_of(set[[2]])
  cat(sprintf(paste0("variance components, %s: %d; abs z median %s, ",
                     "max %s; > 1: %d\n"), set[[1]], nrow(S),
              f3(stats::median(abs(S$z))), f3(max(abs(S$z))),
              sum(abs(S$z) > 1)))
}
S <- sds_of(ok2)
S <- S[abs(S$z) > 1, ]
cat("variance components past abs z 1, both converged:\n")
for (i in seq_len(nrow(S))) {
  cat(sprintf("  %-22s %-30s frmtmb %8s  brms %8s (sd %s)  z %s  %s\n",
              S$id[i], S$name[i], f3(S$frmtmb[i]), f3(S$brms[i]),
              f3(S$brms_sd[i]), f3(S$z[i]),
              if (grepl("^sd_residual", S$name[i])) "residual SD" else
                if (S$frmtmb[i] < S$brms_sd[i])
                  "within one posterior SD of zero" else "group-level SD"))
}

cat("\n## detail: every coefficient with abs z > 1, and every sd\n")
for (id in both$id) {
  d <- detail[[id]]
  big <- d$fx[abs(d$fx$z) > 1, ]
  if (nrow(big) || !is.null(d$sd) || length(d$only_f) || length(d$only_b)) {
    cat("\n###", id, "\n")
  }
  if (nrow(big)) {
    for (j in seq_len(nrow(big))) {
      cat(sprintf("  b  %-28s frmtmb %10s  brms %10s (sd %s)  z %s\n",
                  big$name[j], f3(big$frmtmb[j]), f3(big$brms[j]),
                  f3(big$brms_sd[j]), f3(big$z[j])))
    }
  }
  if (!is.null(d$sd)) {
    for (j in seq_len(nrow(d$sd))) {
      cat(sprintf("  sd %-28s frmtmb %10s  brms %10s (sd %s)  z %s\n",
                  d$sd$name[j], f3(d$sd$frmtmb[j]), f3(d$sd$brms[j]),
                  f3(d$sd$brms_sd[j]), f3(d$sd$z[j])))
    }
  }
  if (length(d$only_f)) cat("  names only in frmtmb:",
                            paste(d$only_f, collapse = ", "), "\n")
  if (length(d$only_b)) cat("  names only in brms:  ",
                            paste(d$only_b, collapse = ", "), "\n")
}
cat("\n## models not compared, and why\n\n")
for (id in setdiff(R$id, both$id)) {
  d <- detail[[id]]
  cat(sprintf("- %s: frmtmb %s%s; brms %s%s\n", id, R$frmtmb[R$id == id],
              if (nzchar(d$f_msg %||% ""))
                paste0(" (", substr(d$f_msg, 1, 100), ")") else "",
              R$brms[R$id == id],
              if (nzchar(d$b_msg %||% ""))
                paste0(" (", substr(d$b_msg, 1, 100), ")") else ""))
}
saveRDS(list(table = R, detail = detail),
        file.path(PORT_OUT, "plausibility.rds"))
