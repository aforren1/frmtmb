# Short check 2 of lane wcheck: the sample route's verdict per link.
# A posterior's integrability does not depend on the parameterization:
# with the Jacobian included, int p(x) L(x) dx is the same over eta as
# over x. Near 0 it needs order + 1 - E > 0, toward infinity
# p - 1 - E > 0, for every link. The MAP objective (a density in eta)
# carries the link's own power j. This script prints the check's verdict
# on both routes for each link, next to that condition.
#   Rscript dev/reviews/wcheck-short2-routes.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
cw <- frmtmb:::check_weights_unbounded
verdict <- function(fr, pr, route) {
  tryCatch({cw(fr$spec, fr, pr, exact = TRUE, route = route); "pass"},
           error = function(e) {
             m <- conditionMessage(e)
             if (grepl("improper", m)) "IMPROPER" else
               if (grepl("no maximum", m)) "REFUSED" else
                 paste("ERROR", substr(m, 1, 80))
           })
}
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n))
d$y <- 0.5 + 0.4 * d$x + rnorm(n)
mu <- exp(0.5 + 0.3 * d$x + rnorm(n, 0, 0.3))
d$yg <- rgamma(n, shape = 4, rate = 4 / mu)
d$yf <- mu * (-log(runif(n)))^(-1 / 4)
cat("sigma toward 0, normal(0, 1) on sigma (order 0); integrable iff",
    "1 - E > 0\n")
for (E in c(0.5, 1.5)) {
  d$w <- 1 + E / n
  for (ln in c("log", "identity", "softplus", "squareplus")) {
    fr <- q(frm(bf(y | weights(w) ~ x + (1 | obs)),
                family = brmsfamily("gaussian", link_sigma = ln), data = d,
                dry_run = "frame"))
    pr <- set_prior("normal(0, 1)", class = "sigma")
    cat(sprintf("  E %g link %-10s map %-8s sample %-8s | posterior %s\n", E,
                ln, verdict(fr, pr, "map"), verdict(fr, pr, "sample"),
                if (1 - E > 0) "proper" else "improper"))
  }
}
cat("Gamma shape toward infinity, E 12 (growth 6), student_t(nu) on",
    "shape (p = nu + 1); integrable iff p - 1 - 6 > 0\n")
d$w <- 1 + 12 / n
for (nu in c(5.5, 6.5, 7.5)) {
  for (ln in c("log", "identity", "softplus", "squareplus")) {
    fr <- q(frm(bf(yg | weights(w) ~ x + (1 | obs)),
                family = brmsfamily("Gamma", link = "log", link_shape = ln),
                data = d, dry_run = "frame"))
    pr <- set_prior(sprintf("student_t(%g, 0, 2.5)", nu), class = "shape")
    cat(sprintf("  nu %g link %-10s map %-8s sample %-8s | posterior %s\n",
                nu, ln, verdict(fr, pr, "map"), verdict(fr, pr, "sample"),
                if (nu + 1 - 1 - 6 > 0) "proper" else "improper"))
  }
}
cat("frechet nu (default link logm1), E 8 (growth 8), student_t(nu0)\n")
d$w <- 1 + 8 / n
fr <- q(frm(bf(yf | weights(w) ~ x + (1 | obs)), family = frechet(),
            data = d, dry_run = "frame"))
cat("  nu link:", fr$linpreds[[grep("nu", names(fr$linpreds))[1]]]$link$name,
    "\n")
for (nu0 in c(3, 12)) {
  pr <- set_prior(sprintf("student_t(%g, 0, 2.5)", nu0), class = "nu")
  cat(sprintf("  student_t(%g): map %s sample %s\n", nu0,
              verdict(fr, pr, "map"), verdict(fr, pr, "sample")))
}
cat("DONE\n")
