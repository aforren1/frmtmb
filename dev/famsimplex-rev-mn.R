# Reviewer: multinomial and eam/latent multi-primary arms of the parser
# change.  Rscript famsimplex-rev-mn.R lane|base
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm:", arm, " lib:", find.package("frmtmb"), "\n")
run <- function(lab, expr) {
  cat("\n---", lab, "\n")
  r <- tryCatch(suppressMessages(suppressWarnings(expr)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  if (inherits(r, "frmtmb_fit")) {
    cat(sprintf("logLik %.10f\n", as.numeric(logLik(r))))
    print(round(fixef(r)[, "Estimate"], 5))
  } else cat(r, "\n")
}
set.seed(31)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
E <- cbind(0, 0.3 + 0.5 * d$z, -0.2 + 0.6 * d$x)
P <- exp(E) / rowSums(exp(E))
d$N <- 10L
d$Ym <- t(vapply(seq_len(n), function(i) rmultinom(1, 10, P[i, ])[, 1],
                 numeric(3)))
colnames(d$Ym) <- c("a", "b", "c")
mn <- multinomial(K = 3)
run("multinomial ~ x", frm(bf(Ym | trials(N) ~ x), family = mn, data = d))
run("multinomial ~ x, mu3 ~ z", frm(bf(Ym | trials(N) ~ x, mu3 ~ z),
                                    family = mn, data = d))
run("multinomial ~ x, mu2 ~ z (new)", frm(bf(Ym | trials(N) ~ x, mu2 ~ z),
                                          family = mn, data = d))
run("multinomial ~ x, mu2 = 0", frm(bf(Ym | trials(N) ~ x, mu2 = 0),
                                    family = mn, data = d))
cat("\nprimary dpars:", tryCatch(mn$primary_dpars, error = function(e) "?"),
    "\n")
