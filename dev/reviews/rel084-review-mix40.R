# rel084 release review, consolidator's semantic finding (b): drawdir's
# mix40 cumulative mixtures, seeds 8, 10, 27 (and 3 controls), where
# ordstarts' swap reaches a higher maximum and predict() then says the
# covariance could not be recovered from the Hessian. Is the Hessian
# singular at that (correct) maximum, does the fit say so at fit time, and
# does predict() say it clearly? Prints the fit's own conditions (the
# drawdir case builder muffles them), logLik, code, which SEs are lost and
# why, the Hessian's smallest eigenvalues, and predict()'s conditions.
#   REVLIB=r23|r22 Rscript dev/reviews/rel084-review-mix40.R
source("dev/reviews/rel084-review-common.R")
options(width = 170)
fq <- NULL
source("dev/predraw2-cases.R")
ns <- asNamespace("frmtmb")
for (s in c(8, 10, 27, 1, 2, 3)) {
  d <- gen("mixture_cumulative", "interior", s)
  sp <- spec("mixture_cumulative")
  r <- conds_of(frm(sp[[1]], family = sp[[2]], data = d))
  f <- r$v
  cat(sprintf("== seed %d ll %.6f code %d | fit: %s\n", s, as.numeric(logLik(f)),
              f$opt$convergence, cls_line(r)))
  for (x in r$said) if (x$k == "W" || grepl("swap|degenerate", x$msg))
    cat("   ", x$k, x$cls[[1]], ":", substr(gsub("\\s+", " ", x$msg), 1, 260), "\n")
  sdr <- tryCatch(ns$sdr_of(f), error = function(e) NULL)
  if (!is.null(sdr)) {
    lost <- sdr$se_lost
    if (length(lost)) cat("    se_lost:", paste(names(lost), lost, sep = "=", collapse = " "), "\n")
  }
  V <- tryCatch(suppressWarnings(vcov(f)), error = function(e) NULL)
  cat("    vcov finite diag:", if (is.null(V)) "NULL" else
        paste(sum(is.finite(diag(V))), "of", nrow(V)), "\n")
  H <- tryCatch(f$obj$he(f$opt$par), error = function(e) NULL)
  if (!is.null(H)) {
    ev <- sort(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values)
    cat("    Hessian eigen (smallest 3, largest):", paste(signif(ev[1:3], 4), collapse = " "),
        signif(max(ev), 4), "\n")
    print(signif(f$opt$par, 5))
  }
  set.seed(5)
  p <- conds_of(predict(f, ndraws = 200, summary = FALSE))
  cat("    predict:", cls_line(p), "\n")
  for (x in p$said) cat("     ", x$cls[[1]], ":", substr(gsub("\\s+", " ", x$msg), 1, 200), "\n")
}
cat("DONE\n")
