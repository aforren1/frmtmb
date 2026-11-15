# Short check of lane priors2's fix round, item 2 (F1): is the refusal
# after the fit a proof? Case H: poisson("identity") rate per unit of a
# large exposure, gamma(0.5, 1) or beta(0.5, 2) on the rate with lb = 0.
# The log posterior in b is S log b - b (X + rate) + (a - 1) log b + c,
# S = sum(y), X = sum(x): with S + a - 1 > 0 it goes to -Inf at b = 0
# and its maximum is b* = (S + a - 1) / (X + rate), finite. A large X puts
# b* very near the end while the objective stays finite.
# Rscript priors2-short-f1.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
show <- function(lab, expr) {
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  [warning]", substr(conditionMessage(w), 1, 110), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    cat(lab, ": ERROR ", substr(conditionMessage(r), 1, 160), "\n", sep = "")
  } else {
    cat(lab, ": code ", r$opt$convergence, " objective ",
        format(r$opt$objective, digits = 15), " par ",
        paste(format(r$opt$par, digits = 10), collapse = " "), "\n",
        sep = "")
  }
  invisible(r)
}
# hand log posterior (negated, up to the constant frm() also carries)
nlp <- function(b, y, x, a, rate) {
  -(sum(y * log(b * x)) - b * sum(x) + (a - 1) * log(b) - rate * b)
}
set.seed(5)
n <- 80
base_x <- runif(n, 0.5, 2)
for (scale in c(1e5, 1e7, 1e9)) {
  for (S in c(1L, 3L)) {
    h <- data.frame(x = base_x * scale, y = 0L)
    h$y[seq_len(S)] <- 1L
    X <- sum(h$x)
    bstar <- (S - 0.5) / (X + 1)
    cat(sprintf("H gamma(0.5,1): scale %g, S = %d, b* = %.6g, -logpost(b*) = %.10g, at b*/1e3: %.10g, at 1e-300: %.6g\n",
                scale, S, bstar, nlp(bstar, h$y, h$x, 0.5, 1),
                nlp(bstar / 1e3, h$y, h$x, 0.5, 1),
                nlp(1e-300, h$y, h$x, 0.5, 1)))
    show(sprintf("  H lb = 0     scale %g S %d", scale, S),
         frm(y ~ 0 + x, data = h, family = poisson(link = "identity"),
             prior = set_prior("gamma(0.5, 1)", class = "b", coef = "x",
                               lb = 0)))
    show(sprintf("  H lb = b*/1e4 scale %g S %d", scale, S),
         frm(y ~ 0 + x, data = h, family = poisson(link = "identity"),
             prior = set_prior("gamma(0.5, 1)", class = "b", coef = "x",
                               lb = bstar / 1e4)))
    show(sprintf("  H beta(0.5,2) lb = 0 ub = 1 scale %g S %d", scale, S),
         frm(y ~ 0 + x, data = h, family = poisson(link = "identity"),
             prior = set_prior("beta(0.5, 2)", class = "b", coef = "x",
                               lb = 0, ub = 1)))
  }
}
