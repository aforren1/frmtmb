# Reviewer, fix round 1: priors on a cell-means dispersion, exgaussian
# with one scale held, and the asym_laplace refusal.
#   Rscript dev/fuzzfix-rev-r1-extra.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
run <- function(name, expr) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    m <- conditionMessage(r)
    tag <- if (grepl("has no maximum|is improper", m)) "REFUSED" else "ERROR"
    cat(sprintf("ARM %-34s %s | %s\n", name, tag, substr(m, 1, 120)))
  } else {
    cat(sprintf("ARM %-34s ok logLik %.6g | betad %s | conv %s | warnings %d%s\n",
                name, as.numeric(logLik(r)),
                paste(sprintf("%.4g", r$estimates$betad), collapse = " "),
                r$opt$convergence %||% NA, length(ws),
                if (length(ws)) paste0(" (", substr(ws[1], 1, 60), ")") else ""))
  }
}
set.seed(20261021)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                fac = factor(rep(c("a", "b"), 40)),
                w = runif(n, 0.5, 2), g = factor(rep(1:20, each = 4)))
eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- eta + rnorm(n, 0, 0.8)
d$ye <- eta + rnorm(n, 0, 0.5) + rexp(n, 1 / 0.6)
d$w12 <- 1.2
G <- gaussian()
cat("--- sigma ~ 0 + fac (cell means span the constant), prior on b_sigma\n")
fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 0 + fac)
run("no prior", frm(fo, family = G, data = d))
run("normal(0,1) on b, dpar sigma", frm(fo, family = G, data = d,
    prior = set_prior("normal(0, 1)", class = "b", dpar = "sigma")))
cat("--- exgaussian, (1 | obs), w = 1.2\n")
EX <- exgaussian()
run("both scales free", frm(bf(ye | weights(w12) ~ x + (1 | obs)),
                            family = EX, data = d))
run("beta held at 0.6", frm(bf(ye | weights(w12) ~ x + (1 | obs), beta = 0.6),
                            family = EX, data = d))
run("sigma held at 0.5", frm(bf(ye | weights(w12) ~ x + (1 | obs), sigma = 0.5),
                             family = EX, data = d))
run("unweighted", frm(bf(ye ~ x + (1 | obs)), family = EX, data = d))
cat("--- asym_laplace\n")
AL <- asym_laplace()
run("obs, unweighted, Laplace", frm(bf(y ~ x + (1 | obs)), family = AL, data = d))
run("obs, unweighted, quadrature", frm(bf(y ~ x + (1 | obs)), family = AL,
                                       data = d, quadrature = TRUE))
run("obs, w = 1.2, quadrature", frm(bf(y | weights(w12) ~ x + (1 | obs)),
                                    family = AL, data = d, quadrature = TRUE))
run("group of 4, unweighted", frm(bf(y ~ x + (1 | g)), family = AL, data = d))
run("group of 4, w = 1.2", frm(bf(y | weights(w12) ~ x + (1 | g)),
                               family = AL, data = d))
d2 <- d; d2$g2 <- factor(c(1:40, rep(41:60, each = 2)))
run("40 singletons + 20 pairs", frm(bf(y ~ x + (1 | g2)), family = AL, data = d2))
run("obs, sigma ~ 0 + fac", frm(bf(y ~ x + (1 | obs), sigma ~ 0 + fac),
                                family = AL, data = d))
cat("DONE\n")
