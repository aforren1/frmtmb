# Release reviewer of 0.71.0: does fit_outer_hessian()'s kept result
# (lane illse) ever serve a pulled Hessian (lane ci701) to an unpulled
# call, or the reverse? test-review-v29.R's boundary-correlation fit
# from theta = (0, 0, 500), and an rs20-like fit with par_units, each
# call compared with a fresh computation (kept result cleared).
#   Rscript dev/relrev071-cache.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
set.seed(1)
ng <- 12
u <- stats::rnorm(ng, 0, 0.8)
per <- 5
dd <- data.frame(id = factor(rep(paste0("a", seq_len(ng)), each = per)),
                 x = stats::rnorm(ng * per))
dd$y1 <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.4)
dd$y2 <- -1 + 0.2 * dd$x + 0.6 * u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.3)
f <- mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
          bf(y2 ~ x + (1 | q | id)) + gaussian())
fit <- suppressMessages(frm(f, data = dd, start = list(theta = c(0, 0, 500))))
check <- function(fit, lbl) {
  np <- length(fit$opt$par)
  ec <- ns$se_edge_cor(fit, np)
  cat(lbl, ": units", if (is.null(fit$par_units)) "NULL" else
        paste(signif(range(fit$par_units), 3), collapse = ".."),
      "| edge cor", which(ec), "\n")
  env <- fit$obj$env
  fresh <- function(pull) {
    rm(list = intersect(".frmtmb_outer_h", ls(env, all.names = TRUE)),
       envir = env)
    ns$fit_outer_hessian(fit, pull = pull)
  }
  hu <- fresh(NULL)
  hp <- fresh(ec)
  cat("  pulled differs from unpulled:", !identical(hu$H, hp$H),
      " max |dH|", signif(max(abs(hu$H - hp$H)), 3), "\n")
  # sequences of calls on the kept result, each against the fresh one
  seqs <- list(c("u", "p", "u", "p"), c("p", "u", "p", "u"),
               c("p", "p", "u", "u"), c("u", "at", "p", "at_kept"))
  rm(list = intersect(".frmtmb_outer_h", ls(env, all.names = TRUE)),
     envir = env)
  for (s in seqs) {
    out <- character()
    for (k in s) {
      ok <- switch(k,
        u = identical(ns$fit_outer_hessian(fit)$H, hu$H),
        p = identical(ns$fit_outer_hessian(fit, pull = ec)$H, hp$H),
        at = identical(ns$fit_outer_hessian_at(fit, fit$opt$par), hu$H),
        at_kept = {
          r <- ns$fit_outer_hessian_at(fit, fit$opt$par, kept_only = TRUE)
          is.null(r) || identical(r, hu$H)
        })
      out <- c(out, paste0(k, "=", ok))
    }
    cat("  sequence", paste(out, collapse = " "), "\n")
  }
  # the kept key after a pulled call names the pulled point
  invisible(ns$fit_outer_hessian(fit, pull = ec))
  k <- env$.frmtmb_outer_h$key$at
  cat("  key after pulled call equals pulled point:",
      identical(k, unname(as.numeric(ns$se_pull_back(
        if (is.null(fit$par_units) || all(fit$par_units == 1))
          ns$sdr_outer_point(fit$obj) else fit$opt$par, ec)))), "\n")
}
check(fit, "v29 t0=500")
# with response units (scaled branch)
dd2 <- dd
dd2$y1 <- dd2$y1 * 1e3
dd2$y2 <- dd2$y2 * 1e-2
fit2 <- suppressMessages(frm(f, data = dd2,
                             start = list(theta = c(0, 0, 500))))
check(fit2, "v29 rescaled t0=500")
