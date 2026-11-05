# Item 7: rerun two vignette figures and the NEWS figures, from what a READER of
# the page would run, and check the page's spelling reproduces them.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
sys.source(system.file("rl", "rw-delta.R", package = "frmtmb", mustWork = TRUE),
           envir = globalenv())
sq <- function(e) suppressWarnings(suppressMessages(e))
# (1) Exams_2, against my own integrate() on my own density
d <- data.frame(k = c(21, 17, 21, 18, 22, 31, 31, 34, 34, 35, 35, 36, 39, 36, 35),
                n = 40L, zero = 0, id = factor(1:15))
fo <- bf(k | trials(n) ~ 1 + (1 | id), mu2 ~ 0 + offset(zero),
         theta1 ~ 0 + offset(zero))
lap <- sq(frm(fo, family = mixture(binomial, binomial), data = d))
imp <- sq(frm(fo, family = mixture(binomial, binomial), data = d, importance = 512))
ex_nll <- function(p) {
  s <- exp(p[2])
  -sum(vapply(d$k, function(k) {
    f <- function(u) {
      exp(log(0.5) + dbinom(k, 40, plogis(p[1] + u), log = TRUE)) *
        dnorm(u, 0, s) + 0.5 * dbinom(k, 40, 0.5) * dnorm(u, 0, s)
    }
    log(integrate(f, -12 * s, 12 * s, subdivisions = 2000L, rel.tol = 1e-11)$value)
  }, 0))
}
pl <- unname(lap$opt$par)
pi <- unname(imp$opt$par)
cat("Exams_2 Laplace", format(as.numeric(logLik(lap)), digits = 8), "corrected",
    format(as.numeric(logLik(imp)), digits = 8), "mcse",
    format(imp$importance$mcse, digits = 3), "\n")
cat("exact at its estimates", format(-ex_nll(pi), digits = 8), " at Laplace's",
    format(-ex_nll(pl), digits = 8), "\n")
ref <- nlminb(pl, ex_nll)
cat("exact max", format(-ref$objective, digits = 8), "at", format(ref$par, digits = 6),
    " Laplace at", format(pl, digits = 6), " corrected at", format(pi, digits = 6), "\n")
# (2) the page's spelling of dataset 4206, first in a state a reader might be in,
# then with the seed the old page set
form <- bf(choice | reward(pay1, pay2) ~ condition + (1 | p | id),
           beta ~ 1 + (1 | p | id))
fam <- rw_delta(subject = id, trial = trial)
lsd <- function(f) log(VarCorr(f)$id$sd[1, "Estimate"])
set.seed(99)
a <- rl_simulate(rl_bandit_design(40, 20), seed = 4206)[[1L]]
set.seed(4206)
b <- rl_simulate(rl_bandit_design(40, 20), seed = 4206)[[1L]]
cat("page's spelling same data at another session seed:", identical(a, b), "\n")
for (nm in c("page", "seeded")) {
  dd <- if (nm == "page") a else b
  l <- sq(frm(form, family = fam, data = dd))
  cat(nm, "Laplace log sd", format(lsd(l), digits = 5), "\n")
  if (nm == "page" && identical(a, b)) next
  for (N in c(100L, 400L)) {
    f <- sq(frm(form, family = fam, data = dd, importance = N))
    cat(nm, "N", N, "log sd", format(lsd(f), digits = 5), "sd",
        format(exp(lsd(f)), digits = 3), "mcse",
        format(f$importance$mcse, digits = 3), "ess",
        format(f$importance$ess_min, digits = 3), "rounds", f$importance$rounds, "\n")
  }
}
# draws grow: is -0.6501 at 400 settled? the page says what the draws buy
dd <- b
for (N in c(1600L, 6400L)) {
  f <- sq(frm(form, family = fam, data = dd, importance = N))
  cat("N", N, "log sd", format(lsd(f), digits = 5), "mcse",
      format(f$importance$mcse, digits = 3), "ess",
      format(f$importance$ess_min, digits = 3), "\n")
}
# (3) NEWS's p124 figures through the survey's own data
E <- new.env(parent = asNamespace("frmtmb"))
sys.source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/equiv-specs.R", envir = E)
db <- E$eq_dat("binomial", seed = 124)
f <- bf(y | trials(ntr) ~ x + (1 | g))
i <- sq(frm(f, family = binomial(), data = db, importance = 200L))
q <- sq(frm(f, family = binomial(), data = db, quadrature = TRUE))
la <- sq(frm(f, family = binomial(), data = db))
cat("NEWS p124 imp b0", format(i$opt$par[[1]], digits = 7), "ll",
    format(as.numeric(logLik(i)), digits = 10), "mcse",
    format(i$importance$mcse, digits = 3), "rounds", i$importance$rounds, "\n")
cat("quad b0", format(q$opt$par[[1]], digits = 7), "ll",
    format(as.numeric(logLik(q)), digits = 10), " Laplace ll",
    format(as.numeric(logLik(la)), digits = 10), "\n")
# (4) NEWS's gaussian claim, 1e-9
set.seed(1003L)
g <- factor(rep(1:30, each = 3L))
x <- round(rnorm(90L), 3)
a1 <- rnorm(30L)
b1 <- 0.3 * a1 + sqrt(1 - 0.3^2) * rnorm(30L)
dg <- data.frame(y = 0.3 * x + 0.8 * a1[g] + 0.4 * b1[g] * x + 0.2 +
                   rnorm(90L, 0, 0.5), x = x, g = g)
lg <- sq(frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = dg))
for (k in c(4L, 16L, 64L, 200L)) {
  fg <- sq(frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = dg, importance = k))
  cat("gaussian N", k, "ll - Laplace", format(as.numeric(logLik(fg)) - as.numeric(logLik(lg)),
                                              digits = 3),
      "max par diff", format(max(abs(fg$opt$par - lg$opt$par)), digits = 3), "\n")
}
# (5) the profile figure in NEWS: 40 groups of 4, 500 draws
set.seed(11)
ng <- 40L
per <- 4L
g2 <- factor(rep(seq_len(ng), each = per))
