# Item 3631: the joint-density translator (stan_pars_from_fit(),
# tests/testthat/helper-brms.R) on multinomial() with (1 | g), where
# check C fails, against (1 | ID | g), where it passes. Seed 65,
# nbdm_dm_data() of test-nb2dm-brms.R.
#   Rscript dev/hygiene-3631.R <mine|base> > dev/hygiene-out/3631.log
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (identical(arm, "mine")) "C:/Users/adf44/source/r/wt-hygiene-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-hygiene"
Sys.setenv(FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true",
           FRMTMB_STAN_CACHE = file.path(WT, "dev", "stan-cache"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER"))) Sys.setenv(R_MAKEVARS_USER = mk)
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
H <- new.env(parent = asNamespace("frmtmb"))
hp <- commandArgs(TRUE)[2]
if (is.na(hp)) hp <- file.path(WT, "tests", "testthat", "helper-brms.R")
for (h in hp) {
  sys.source(h, envir = H)
}
options(frmtmb.brms_lp_report = TRUE)
set.seed(65)
n <- 300
dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:15, length.out = n)))
E <- cbind(0, 0.3 + 0.6 * dd$x + stats::rnorm(15, 0, 0.4)[dd$g],
           -0.4 - 0.5 * dd$x + stats::rnorm(15, 0, 0.3)[dd$g])
P <- exp(E) / rowSums(exp(E))
dd$N <- sample(4:25, n, TRUE)
Yc <- t(vapply(seq_len(n), function(i) {
  as.numeric(stats::rmultinom(1, dd$N[i], P[i, ]))
}, numeric(3)))
colnames(Yc) <- c("y1", "y2", "y3")
dd$Yc <- Yc
show <- function(label, f_frm, f_brms) {
  cat("\n== ", label, "\n", sep = "")
  fit <- frm(f_frm, family = multinomial(), data = dd)
  re <- ranef(fit)
  cat("ranef blocks:", paste(names(re), collapse = " | "), "\n")
  for (nm in names(re)) cat("  ", nm, ": cols", colnames(re[[nm]]), "\n")
  vc <- unclass(frmtmb:::varcorr_matrices(fit))
  cat("VarCorr blocks:", paste(names(vc), collapse = " | "), "\n")
  for (nm in names(vc)) cat("  ", nm, ": cols", colnames(vc[[nm]]), "\n")
  prior <- H$brms_flat_prior(f_brms, data = dd, family = brms::multinomial())
  rtab <- H$brms_ranef_table(f_brms, dd, brms::multinomial(), prior)
  print(rtab[, intersect(c("id", "group", "gn", "coef", "cn", "resp",
                           "dpar", "nlpar", "cor"), names(rtab))])
  r <- with_reporter("check", H$brms_lp_check(f_brms, brms::multinomial(), dd,
                                              fit, joint = TRUE))
  cat(sprintf("RESULT %s const %.4g grad %.3g\n", label, r$measured_const,
              r$max_grad))
  invisible(r)
}
r1 <- tryCatch(show("(1 | g)", Yc | trials(N) ~ x + (1 | g),
                    brms::bf(Yc | trials(N) ~ x + (1 | g))),
               error = function(e) cat("ERROR", conditionMessage(e), "\n"))
r2 <- tryCatch(show("(1 | ID | g)", Yc | trials(N) ~ x + (1 | ID | g),
                    brms::bf(Yc | trials(N) ~ x + (1 | ID | g))),
               error = function(e) cat("ERROR", conditionMessage(e), "\n"))
