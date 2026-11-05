# The latent residual beside other model features: L1 against brms 2.23.0's
# log_prob, through the survey's instrument, as dev/latres-l1.R does, and a
# record of whatever either side refuses. Seed 167.
#
# Usage: Rscript dev/latres-l1x.R <worker> <workers> <outfile>
args <- commandArgs(TRUE)
worker <- as.integer(args[[1]])
workers <- as.integer(args[[2]])
outf <- args[[3]]
setwd("C:/Users/adf44/source/r/frmtmb-wt-latres")
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(rstan)
})
# The instrument reads the objective, which frm_sample() reads too, on the forms
# frm() refuses for its Laplace value.
assignInNamespace("latres_refuse_laplace", function(spec, frame) NULL, "frmtmb")
assignInNamespace("latres_refuse_margin", function(fit) NULL, "frmtmb")
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- function(a, b) if (is.null(a)) b else a
sys.source("dev/equiv-lib.R", envir = hel)

set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:30))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
d$z <- rnorm(n)
d$h <- factor(rep(1:6, length.out = n))
e <- as.vector(replicate(30, arima.sim(list(ar = 0.5), 5, sd = 0.6)))
eta <- 0.8 + 0.4 * d$x + e
d$count <- rpois(n, exp(eta))
d$count2 <- rpois(n, exp(0.4 + 0.2 * d$x + rev(e)))
d$bin <- rbinom(n, 1, plogis(eta - 0.8))
d$pos <- rgamma(n, 4, 4 / exp(eta))
d$ord <- factor(findInterval(eta + rlogis(n), c(0, 1, 2)) + 1L, ordered = TRUE)
d$w <- runif(n, 0.2, 0.9)
d$cen <- rbinom(n, 1, 0.1)
d$mo3 <- sample(1:3, n, TRUE)
d$sx <- 0.2
drag <- d[-c(3, 4, 17, 40, 41, 42, 77, 90), ]
drag$ragy <- drag$count

ar <- "ar(t, g, cov = TRUE)"
specs <- list(
  re = list("count ~ x + (1 | h) + ar(t, g, cov = TRUE)", "poisson()"),
  reg = list("count ~ x + (1 | g) + ar(t, g)", "poisson()"),
  weights = list("count | weights(w) ~ x + ar(t, g, cov = TRUE)", "poisson()"),
  cens = list("count | cens(cen) ~ x + ar(t, g, cov = TRUE)", "poisson()"),
  trunc = list("count | trunc(lb = 0) ~ x + ar(t, g, cov = TRUE)", "poisson()"),
  censf = list("count | cens(cen) ~ x + ar(t, g)", "poisson()"),
  rag = list("count ~ x + ar(t, g, cov = TRUE)", "poisson()", drag),
  ragma = list("count ~ x + ma(t, g, cov = TRUE)", "poisson()", drag),
  ragarma = list("count ~ x + arma(t, g, cov = TRUE)", "poisson()", drag),
  ragcosy = list("count ~ x + cosy(t, g)", "poisson()", drag),
  ragunstr = list("count ~ x + unstr(t, g)", "poisson()", drag),
  ragar2 = list("count ~ x + ar(t, g, p = 2)", "poisson()", drag),
  nogroup = list("count ~ x + ar(tt, cov = TRUE)", "poisson()",
                 transform(d[seq_len(40L), ], tt = seq_len(40L))),
  notime = list("count ~ x + ar(gr = g, cov = TRUE)", "poisson()"),
  sig = list("bf(pos ~ x + ar(t, g, cov = TRUE), sigma ~ z)", "lognormal()"),
  shape = list("bf(pos ~ x + ar(t, g), shape ~ z)", "Gamma(link = 'log')"),
  zi = list("bf(count ~ x + ar(t, g, cov = TRUE), zi ~ z)",
            "zero_inflated_poisson()"),
  smooth = list("count ~ s(z, k = 5) + ar(t, g, cov = TRUE)", "poisson()"),
  mo = list("count ~ mo(mo3) + ar(t, g, cov = TRUE)", "poisson()"),
  cs = list("ord ~ x + cs(z) + ar(t, g, cov = TRUE)", "sratio()"),
  acat = list("ord ~ x + ar(t, g)", "acat()"),
  mv = list("bf(count ~ x + ar(t, g, cov = TRUE)) + bf(count2 ~ x + ma(t, g, cov = TRUE))",
            "poisson()"),
  probit = list("bin ~ x + ar(t, g)", "bernoulli('probit')"),
  psqrt = list("count ~ x + ar(t, g)", "poisson('sqrt')"),
  ident = list("count ~ x + ar(t, g, cov = TRUE)", "poisson('identity')"),
  autocorarg = list("bf(count ~ x, autocor = ~ ar(t, g, cov = TRUE))",
                    "poisson()"),
  acf = list("bf(count ~ x) + acformula(~ ar(t, g, cov = TRUE))", "poisson()"))
if (length(args) > 3L) specs <- specs[strsplit(args[[4]], ",")[[1]]]
idx <- seq_along(specs)
idx <- idx[(idx - 1L) %% workers == worker - 1L]
if (!file.exists(outf)) {
  cat(paste(c("spec", "const", "const_up", "const_dn", "grad", "loglik",
              "conv", "note"), collapse = "\t"), "\n", sep = "", file = outf)
}
mkf <- function(txt, pkg) {
  e <- new.env(parent = asNamespace(pkg))
  if (grepl("^bf[(]", txt)) return(eval(parse(text = txt), e))
  eval(parse(text = paste0("bf(", txt, ")")), e)
}
for (i in idx) {
  s <- specs[[i]]
  nm <- names(specs)[i]
  dd <- if (length(s) > 2L) s[[3]] else d
  rec <- rep("", 7)
  res <- tryCatch({
    fam_f <- hel$eq_eval(s[[2]], "frmtmb")
    fam_b <- hel$eq_eval(s[[2]], "brms")
    fit <- suppressWarnings(frm(mkf(s[[1]], "frmtmb"), family = fam_f,
                                data = dd))
    bform <- mkf(s[[1]], "brms")
    prior <- hel$brms_flat_prior(bform, data = dd, family = fam_b)
    refit <- function(k) {
      suppressWarnings(frm(mkf(s[[1]], "frmtmb"), family = fam_f, data = dd,
        control = frmtmb_control(optimizer = hel$eq_step_optimizer(c(1, -1)[k]),
                                 restarts = 0, check_se = "ignore")))
    }
    l1 <- hel$eq_l1(list(data = dd), bform, fam_b, dd, fit, prior,
                    joint = TRUE, refit = refit)
    ex <- function(k) {
      if (length(l1$extra) >= k && !is.null(l1$extra[[k]]$const)) {
        sprintf("%.6g", l1$extra[[k]]$const)
      } else "NA"
    }
    c(sprintf("%.6g", l1$at$const), ex(1L), ex(2L),
      sprintf("%.3g", l1$at$grad), sprintf("%.6f", logLik(fit)),
      fit$opt$convergence,
      paste(vapply(l1$extra, function(e) substr(e$error %||% "", 1, 60), ""),
            collapse = ";"))
  }, error = function(e) {
    c(rep("", 6), gsub("[\t\r\n]+", " ", substr(conditionMessage(e), 1, 300)))
  })
  cat(paste(c(nm, res), collapse = "\t"), "\n", sep = "", file = outf,
      append = TRUE)
}
