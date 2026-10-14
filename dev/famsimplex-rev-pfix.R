# Reviewer: a fixed value on a location dpar of a multi-location family.
#   Rscript famsimplex-rev-pfix.R lane|base   -> dev/famsimplex-rev-pfix-<arm>.log
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm:", arm, " lib:", find.package("frmtmb"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-famsimplex/dev/famsimplex-sim.R")
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
d$f <- factor(c("a", "b", "c")[apply(P, 1, function(p) sample(3, 1, prob = p))])
cl <- rbinom(n, 1, 0.4)
d$y <- ifelse(cl == 1, 2 + 1 * d$x, -1 + 0.5 * d$z) + rnorm(n, 0, 0.7)
mix <- mixture(gaussian(), gaussian())
run("mixture y ~ x", frm(bf(y ~ x), family = mix, data = d))
run("mixture y ~ x, mu2 = 0", frm(bf(y ~ x, mu2 = 0), family = mix, data = d))
run("mixture y ~ x, mu1 = 0", frm(bf(y ~ x, mu1 = 0), family = mix, data = d))
run("mixture y ~ x, mu2 = \"mu1\"",
    frm(bf(y ~ x, mu2 = "mu1"), family = mix, data = d))
run("categorical f ~ x", frm(bf(f ~ x), family = categorical(), data = d))
run("categorical f ~ x, muc = 0",
    frm(bf(f ~ x, muc = 0), family = categorical(), data = d))
run("categorical f ~ x, mub = 0",
    frm(bf(f ~ x, mub = 0), family = categorical(), data = d))
if (arm == "lane") {
  dd <- fs_data_dirichlet(seed = 11)
  run("dirichlet Y ~ x", frm(bf(Y ~ x), family = dirichlet(), data = dd))
  run("dirichlet Y ~ x, muc = 0",
      frm(bf(Y ~ x, muc = 0), family = dirichlet(), data = dd))
  run("dirichlet Y ~ x, mub = 0",
      frm(bf(Y ~ x, mub = 0), family = dirichlet(), data = dd))
  run("dirichlet Y ~ x, phi = 5",
      frm(bf(Y ~ x, phi = 5), family = dirichlet(), data = dd))
  dl <- fs_data_logistic_normal(seed = 13)
  run("logistic_normal Y ~ x, sigmac = 1",
      frm(bf(Y ~ x, sigmac = 1), family = logistic_normal(), data = dl))
  run("logistic_normal Y ~ x, muc = 0",
      frm(bf(Y ~ x, muc = 0), family = logistic_normal(), data = dl))
}
