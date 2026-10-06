# Reviewer, re-check: defect 9. frm_allfit() on a nonlinear fit that
# needs start =. Own construction: brms_nonlinear's fit_loss (the
# vignette model with the lane's start patch), then frm_allfit(). The
# control case is a linear fit, where start = NULL is harmless.
#
#   Rscript dev/vigport-rev2-allfit.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "\n")
loss <- brms::loss
fit_loss <- frm(bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
                   ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE),
                data = loss, family = gaussian(),
                start = list(beta = c(5000, 1, 45)))
cat("fit_loss: code", fit_loss$opt$convergence, " logLik",
    format(as.numeric(logLik(fit_loss)), digits = 10), "\n")
show_allfit <- function(f, tag) {
  a <- tryCatch(frm_allfit(f), error = function(e) {
    cat(tag, "frm_allfit ERROR:", conditionMessage(e), "\n"); NULL
  })
  if (is.null(a)) return(invisible())
  cat("\n==", tag, "\n")
  for (nm in names(a$fits)) {
    g <- a$fits[[nm]]
    cat(sprintf("  %-13s %s\n", nm, if (is.null(g)) "NULL (refit failed)"
                else sprintf("logLik %.4f  code %s", as.numeric(logLik(g)),
                             g$opt$convergence)))
  }
}
show_allfit(fit_loss, "nonlinear fit_loss")
# control: a linear mixed model, where the default start is fine
d <- brms::epilepsy
fl <- frm(count ~ zAge + zBase * Trt + (1 | patient), data = d,
          family = poisson())
show_allfit(fl, "linear control (epilepsy poisson)")
