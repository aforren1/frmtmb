# Does one ar() parameter take two priors, one per spelling? The
# resolved entries, and the objective's penalty against a hand sum.
# Arg 1: "before" or "after". Log: dev/acprior-out/double-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-acprior-lib",
         "C:/Users/adf44/source/r/rellib-r13",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
set.seed(3)
n <- 60
d <- data.frame(x = rnorm(n), t = rep(1:20, 3), g = rep(1:3, each = 20))
e <- numeric(n)
for (i in seq_len(n)) {
  e[i] <- rnorm(1) + if (d$t[i] > 1) 0.6 * e[i - 1] else 0
}
d$y <- 1 + 0.5 * d$x + e
try_it <- function(lab, f, pr) {
  cat("\n---", lab, "\n")
  r <- tryCatch({
    fit <- frm(f, data = d, prior = pr)
    ri <- frmtmb:::resolve_prior_input(fit, pr)
    for (en in ri$entries) {
      cat("  entry", en$comp, paste(en$idx, collapse = ","), en$dist$kind,
          if (identical(en$dist$kind, "trans")) en$dist$inner$kind, "\n")
    }
    cat("  thetaac at fit:", format(fit$estimates$thetaac, digits = 10),
        "\n")
  }, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
}
for (cv in c(FALSE, TRUE)) {
  f <- y ~ x + ar(time = t, gr = g, cov = cv)
  try_it(paste("cov =", cv, ": ar normal(0.5, 0.2) + theta normal(0, 1)"),
         f, c(set_prior("normal(0.5, 0.2)", class = "ar"),
              set_prior("normal(0, 1)", class = "theta",
                        coef = "thetaac_1")))
  try_it(paste("cov =", cv, ": theta normal(0, 1) + ar normal(0.5, 0.2)"),
         f, c(set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1"),
              set_prior("normal(0.5, 0.2)", class = "ar")))
  try_it(paste("cov =", cv, ": ar normal(0.5, 0.2) only"), f,
         set_prior("normal(0.5, 0.2)", class = "ar"))
  try_it(paste("cov =", cv, ": theta bounds only"), f,
         set_prior("", class = "theta", coef = "thetaac_1", lb = -0.5))
}
try_it("cosy: cosy beta(2, 2) + theta normal(0, 1)",
       y ~ x + cosy(time = t, gr = g),
       c(set_prior("beta(2, 2)", class = "cosy"),
         set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1")))
try_it("unstr: cortime lkj(2) + theta normal(0, 1)",
       y ~ x + unstr(time = t, gr = g),
       c(set_prior("lkj(2)", class = "cortime"),
         set_prior("normal(0, 1)", class = "theta", coef = "thetaac_1")))
