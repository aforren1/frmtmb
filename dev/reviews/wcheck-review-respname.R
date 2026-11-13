# Review of lane wcheck: does the per-response prior lookup match a
# prior's resp to a response whose name brms rewrites (y_a -> ya)?
# A bounded model (lognormal(0, 1) on each response's sigma) must not be
# refused; the same model on y / y2 is not (dev/wcheck-cases.R).
#   Rscript dev/reviews/wcheck-review-respname.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
set.seed(20261021)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n), w = runif(n, 0.5, 2))
d$y_a <- 0.5 + 0.4 * d$x + rnorm(n)
d$y.b <- d$y_a + rnorm(n)
d$ya2 <- d$y_a
d$yb2 <- d$y.b
P <- function(...) set_prior(...)
out <- function(lab, expr) {
  r <- tryCatch({f <- suppressWarnings(suppressMessages(expr))
    sprintf("fit logLik %.6g betad %s conv %d", as.numeric(logLik(f)),
            paste(format(f$estimates$betad, digits = 4), collapse = " "),
            f$opt$convergence)},
    error = function(e) substr(gsub("\\s+", " ", conditionMessage(e)), 1, 200))
  cat(sprintf("%-40s %s\n", lab, r))
}
fo1 <- mvbf(bf(y_a | weights(w) ~ x + (1 | obs)),
            bf(y.b | weights(w) ~ x + (1 | obs)))
fo2 <- mvbf(bf(ya2 | weights(w) ~ x + (1 | obs)),
            bf(yb2 | weights(w) ~ x + (1 | obs)))
out("names ya2 / yb2, resp ya2 / yb2", frm(fo2, family = gaussian(), data = d,
    prior = c(P("lognormal(0, 1)", class = "sigma", resp = "ya2"),
              P("lognormal(0, 1)", class = "sigma", resp = "yb2"))))
out("names y_a / y.b, resp ya / yb", frm(fo1, family = gaussian(), data = d,
    prior = c(P("lognormal(0, 1)", class = "sigma", resp = "ya"),
              P("lognormal(0, 1)", class = "sigma", resp = "yb"))))
out("names y_a / y.b, resp y_a / y.b", frm(fo1, family = gaussian(), data = d,
    prior = c(P("lognormal(0, 1)", class = "sigma", resp = "y_a"),
              P("lognormal(0, 1)", class = "sigma", resp = "y.b"))))
cat("DONE\n")
