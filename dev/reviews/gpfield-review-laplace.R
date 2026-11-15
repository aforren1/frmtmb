# Review of lane gpfield, frmtmb_gp_laplace: false alarms on well-behaved
# exact gp() fits, and misses on white-field bernoulli fits.
# Usage: Rscript gpfield-review-laplace.R <set>   set: ok | white
# ok: 20 fits (poisson, binomial, bernoulli, negbinomial, Gamma, student
# with a smooth field), each warned row read by importance sampling.
# white: 40 bernoulli y ~ gp(x) seeds, small n; every fit whose length
# scale ends under half the median gap is read by importance sampling, and
# whether it warned.
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
set <- commandArgs(TRUE)[1]
gf_libs("lane")
ns <- asNamespace("frmtmb")
tkof <- function(fit) tryCatch(ns$gp_laplace_tk(fit), error = function(e) NA)
row_out <- function(tag, d, fml, fam) {
  t0 <- proc.time()[["elapsed"]]
  q <- gf_fit(bf(fml), data = d, family = fam)
  w <- proc.time()[["elapsed"]] - t0
  f <- q$fit
  if (inherits(f, "error")) {
    cat("LAP", tag, "ERROR", conditionMessage(f), "\n"); return(invisible())
  }
  bk <- Filter(function(b) b$covstruct == "gp", f$frame$re_blocks)[[1]]
  ls <- exp(f$estimates$theta[bk$theta_idx[2]])
  mn <- gf_mednn(d$x)
  lap <- any(grepl("Laplace approximation", q$warn))
  tk <- tkof(f)
  is <- ""
  if (lap || ls < 0.5 * mn || set == "ok") {
    r <- tryCatch(gf_is(f$obj, f$opt$par, ndraw = 10000L),
                  error = function(e) NULL)
    if (!is.null(r)) {
      is <- sprintf("is-lap %.4f (halves %.4f %.4f) ess %.0f ev_min %.3g",
                    r[["is"]] - r[["laplace"]], r[["h1"]] - r[["laplace"]],
                    r[["h2"]] - r[["laplace"]], r[["ess"]], r[["ev_min"]])
    }
  }
  cat(sprintf("LAP %-26s ll %.4f conv %d wall %.1f | ls/mednn %.3g | tk %.4f | lapwarn %d | %s | warn [%s]\n",
              tag, as.numeric(logLik(f)), f$opt$convergence, w, ls / mn, tk,
              lap, is, paste(unique(sub(":.*", "", q$warn)),
                             collapse = ";")))
}
if (set == "ok") {
  i <- 0L
  for (n in c(50, 100, 150)) for (cv in c("sin", "wig")) {
    i <- i + 1L; set.seed(9700 + i)
    x <- sort(round(runif(n, 0, 10), 3))
    f <- if (cv == "sin") sin(x) else sin(2 * x) + 0.2 * x
    row_out(sprintf("pois %s n%d", cv, n),
            data.frame(x, y = rpois(n, exp(1.8 + 0.6 * f))), y ~ gp(x),
            poisson())
  }
  for (n in c(60, 120)) for (cv in c("sin", "wig")) {
    i <- i + 1L; set.seed(9700 + i)
    x <- sort(round(runif(n, 0, 10), 3))
    f <- if (cv == "sin") sin(x) else sin(2 * x) + 0.2 * x - 1
    row_out(sprintf("binom10 %s n%d", cv, n),
            data.frame(x, tr = 10, y = rbinom(n, 10, plogis(f))),
            y | trials(tr) ~ gp(x), binomial())
  }
  for (n in c(150, 300)) for (cv in c("sin", "wig")) {
    i <- i + 1L; set.seed(9700 + i)
    x <- sort(round(runif(n, 0, 10), 3))
    f <- if (cv == "sin") 1.5 * sin(x) else sin(2 * x) + 0.2 * x - 1
    row_out(sprintf("bern %s n%d", cv, n),
            data.frame(x, y = rbinom(n, 1, plogis(f))), y ~ gp(x),
            bernoulli())
  }
  for (n in c(60, 120)) {
    i <- i + 1L; set.seed(9700 + i)
    x <- sort(round(runif(n, 0, 10), 3))
    row_out(sprintf("negbin sin n%d", n),
            data.frame(x, y = rnbinom(n, size = 5, mu = exp(2 + sin(x)))),
            y ~ gp(x), negbinomial())
  }
  for (n in c(60, 120)) {
    i <- i + 1L; set.seed(9700 + i)
    x <- sort(round(runif(n, 0, 10), 3))
    row_out(sprintf("gamma sin n%d", n),
            data.frame(x, y = rgamma(n, 5, 5 / exp(sin(x)))),
            y ~ gp(x), Gamma(link = "log"))
  }
  for (n in c(60, 120)) {
    i <- i + 1L; set.seed(9700 + i)
    x <- sort(round(runif(n, 0, 10), 3))
    row_out(sprintf("student sin n%d", n),
            data.frame(x, y = sin(x) + 0.3 * rt(n, 4)), y ~ gp(x),
            student())
  }
} else {
  for (s in 1:40) {
    set.seed(9800 + s)
    n <- c(30, 50, 80, 120)[1 + s %% 4]
    x <- sort(round(runif(n, 0, 10), 3))
    f <- if (s %% 2) 2 * sin(3 * x) else sin(x)
    row_out(sprintf("bernw s%d n%d", s, n),
            data.frame(x, y = rbinom(n, 1, plogis(1.5 * f))), y ~ gp(x),
            bernoulli())
  }
}
