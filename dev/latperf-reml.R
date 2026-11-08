# Lane latperf: the sparse joint covariance on a REML fit, where `beta` is an
# inner parameter, against the dense route on the same fit. Seed 11, the
# generator of dev/latperf-agree.R, gaussian is no latent residual, so a
# lognormal and a Gamma(log) response.
.libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
src <- readLines("dev/latperf-agree.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^cases <- ", src) - 1L)]))
d <- mk(11, 30, 12)
d$pos <- exp(0.5 + 0.4 * d$x + log(d$cnt + 1) / 3 + stats::rnorm(nrow(d), 0,
                                                                    0.3))
for (fam in list(lognormal(), Gamma(link = "log"))) {
  out <- tryCatch({
    fit <- suppressWarnings(suppressMessages(frm(
      pos ~ x + (1 | s) + ar(t, g, cov = TRUE), family = fam, data = d,
      REML = TRUE)))
    jc <- frmtmb:::get_joint_cov(fit)
    a <- frm_linpred(fit, se.fit = TRUE)$se.fit
    fd <- dense_copy(fit)
    b <- with_dense(frm_linpred(fd, se.fit = TRUE)$se.fit)
    r <- rel(a, b)
    sprintf("sparse %s beta inner %s | linpred %.2e (%d, mis %d)",
            frmtmb:::jc_is_sparse(jc),
            "beta" %in% names(fit$obj$env$par)[fit$obj$env$random],
            r[1], r[2], r[3])
  }, error = function(e) paste("ERROR", conditionMessage(e)))
  cat("RESULT REML", fam$family, out, "\n")
}
