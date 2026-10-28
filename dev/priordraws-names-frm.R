# frm_sample(sample_prior = "yes")'s prior draw names on the models of
# dev/priordraws-brms-names.R, against brms 2.23.0's names read from that
# script's logs (dev/priordraws-log/brms-names-[123].txt). The data are
# built as there. Usage: Rscript dev/priordraws-names-frm.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/wt-priordraws-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), "\n")
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60), g = factor(rep(1:6, 10)),
                h = factor(rep(1:5, each = 12)), m = sample(1:4, 60, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(60)
d$yb <- rbinom(60, 1, 0.5)
d$yp <- rbeta(60, 2, 2)
d$yc <- rpois(60, 3)
d$yo <- pmin(d$yc, 4) + 1
d$t <- rep(1:10, 6)
# brms's names, from the logs
brms_names <- list()
for (f in sprintf("dev/priordraws-log/brms-names-%d.txt", 1:3)) {
  L <- readLines(f)
  cur <- NULL
  for (l in L) {
    if (grepl("^===== ", l)) cur <- sub("^===== (.*) =====$", "\\1", l)
    if (grepl("^prior_draws names:", l)) {
      v <- strsplit(trimws(sub("^prior_draws names:", "", l)), " ")[[1]]
      brms_names[[cur]] <- v[nzchar(v)]
    }
  }
}
P <- function(...) set_prior(...)
cases <- list(
  A = list(bf(y ~ x + (1 | g)), gaussian()),
  B = list(bf(y ~ x + z + (1 + x | g)), gaussian(),
           P("normal(0, 1)", class = "b") +
             P("normal(0, 2)", class = "b", coef = "z") +
             P("normal(1, 3)", class = "Intercept") +
             P("exponential(2)", class = "sigma") +
             P("lkj(2)", class = "cor") +
             P("cauchy(0, 1)", class = "sd", group = "g", coef = "x")),
  B2 = list(bf(y ~ x + (1 + x | g) + (1 | h)), gaussian(),
            P("normal(0, 1)", class = "sd")),
  B3 = list(bf(y ~ x + z), gaussian(),
            P("normal(0, 2)", class = "b", coef = "z")),
  C2 = list(bf(y ~ x + z), gaussian(),
            P("normal(0, 1)", class = "b", lb = 0, ub = 1)),
  H = list(bf(y ~ 0 + Intercept + x), gaussian(),
           P("normal(0, 1)", class = "b")),
  D = list(bf(y ~ mo(m)), gaussian()),
  D2 = list(bf(y ~ mo(m) + x), gaussian(), P("normal(0, 1)", class = "b")),
  I = list(bf(y ~ x, sigma ~ z), gaussian(),
           P("normal(0, 1)", class = "b") +
             P("normal(0, 1)", class = "b", dpar = "sigma")),
  O = list(bf(y ~ x) + bf(z ~ x) + set_rescor(TRUE), NULL,
           P("normal(0, 1)", class = "b", resp = "y")),
  R = list(bf(y ~ x), mixture(gaussian(), gaussian()),
           P("normal(0, 3)", class = "Intercept", dpar = "mu1") +
             P("normal(1, 3)", class = "Intercept", dpar = "mu2")),
  P = list(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE), gaussian(),
           P("normal(1, 2)", nlpar = "a") + P("normal(0, 1)", nlpar = "b")),
  N = list(bf(yo ~ x), cumulative(), P("normal(0, 1)", class = "b")),
  M = list(bf(y ~ s(x)), gaussian()),
  E = list(bf(yp ~ x), Beta()),
  K = list(bf(y ~ x), student()),
  S = list(bf(yc ~ x), hurdle_poisson())
)
res <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  ds <- tryCatch(suppressWarnings(suppressMessages(frm_sample(
    cs[[1]], data = d, family = cs[[2]], prior = if (length(cs) > 2) cs[[3]],
    sample_prior = "yes", chains = 1, iter = 200, seed = 1, refresh = 0))),
    error = function(e) e)
  fn <- if (inherits(ds, "error")) paste("ERROR:", conditionMessage(ds)) else
    names(prior_draws(ds))
  bn <- brms_names[[nm]]
  ok <- identical(fn, bn)
  res[[nm]] <- ok
  cat(sprintf("%-3s %s\n  frmtmb: %s\n  brms:   %s\n", nm,
              if (ok) "SAME" else "DIFF", paste(fn, collapse = " "),
              paste(bn, collapse = " ")))
}
cat("\nSUMMARY same", sum(unlist(res)), "of", length(res), "\n")
