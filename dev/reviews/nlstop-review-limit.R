# Review of lane nlstop, frmtmb_nl_limit: false positives on well
# identified nonlinear fits (each has a proper interior maximum at the
# true values' scale) and false negatives on other limits of a body.
# Usage: Rscript dev/reviews/nlstop-review-limit.R <model,...> <n,...> <seeds>
# Each row: the warnings' classes, the messages' classes, logLik, code,
# and the estimates.
source("dev/nlstop-common.R")
args <- commandArgs(TRUE)
models <- strsplit(args[1], ",")[[1]]
ns_ <- as.integer(strsplit(args[2], ",")[[1]])
seeds <- eval(parse(text = args[3]))
gen <- function(m, n, seed) {
  set.seed(seed)
  switch(m,
    mm = {
      d <- data.frame(x = runif(n, 0, 10))
      d$y <- 10 * d$x / (1 + d$x) + rnorm(n, 0, 0.5)
      list(f = bf(y ~ Vm * x / (K + x), Vm + K ~ 1, nl = TRUE), d = d,
           fam = gaussian(), pr = NULL)
    },
    decay = {
      d <- data.frame(x = runif(n, 0, 8))
      d$y <- 5 * exp(-0.7 * d$x) + rnorm(n, 0, 0.3)
      list(f = bf(y ~ A * exp(-k * x), A + k ~ 1, nl = TRUE), d = d,
           fam = gaussian(), pr = NULL)
    },
    logis = {
      d <- data.frame(x = runif(n, 0, 10))
      d$y <- 10 / (1 + exp((5 - d$x) / 1)) + rnorm(n, 0, 0.5)
      list(f = bf(y ~ Asym / (1 + exp((xmid - x) / scal)),
                  Asym + xmid + scal ~ 1, nl = TRUE), d = d,
           fam = gaussian(), pr = NULL)
    },
    emax = {
      d <- data.frame(x = runif(n, 0, 20))
      d$y <- 1 + 10 * d$x / (2 + d$x) + rnorm(n, 0, 1)
      list(f = bf(y ~ E0 + Emax * x / (EC50 + x), E0 + Emax + EC50 ~ 1,
                  nl = TRUE), d = d, fam = gaussian(),
           pr = set_prior("normal(10, 5)", nlpar = "Emax"))
    },
    nlfs = {
      d <- data.frame(x = runif(n, 0, 5), z = runif(n, -1, 1))
      d$y <- 4 * exp(-0.5 * d$x) +
        rnorm(n, 0, exp(-1 + 0.6 * d$z))
      list(f = bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE) +
             nlf(sigma ~ s0 + s1 * z) + lf(s0 ~ 1, s1 ~ 1), d = d,
           fam = gaussian(), pr = NULL)
    },
    poismm = {
      d <- data.frame(x = runif(n, 0, 10))
      d$y <- rpois(n, 8 * d$x / (1.5 + d$x))
      list(f = bf(y ~ log(Vm * x / (K + x)), Vm + K ~ 1, nl = TRUE),
           d = d, fam = poisson(), pr = NULL)
    },
    bernlogis = {
      d <- data.frame(x = runif(n, 0, 10))
      d$y <- rbinom(n, 1, plogis((d$x - 4) / 1.2))
      list(f = bf(y ~ (x - m) / s, m + s ~ 1, nl = TRUE), d = d,
           fam = bernoulli(), pr = NULL)
    },
    hill = {
      d <- data.frame(x = runif(n, 0.1, 10))
      d$y <- 10 * d$x^2 / (3^2 + d$x^2) + rnorm(n, 0, 0.5)
      list(f = bf(y ~ Emax * x^h / (EC50^h + x^h), Emax + EC50 + h ~ 1,
                  nl = TRUE), d = d, fam = gaussian(), pr = NULL)
    },
    # false negatives: limits the warning is meant for
    fn_explin = {
      d <- data.frame(x = runif(n, 0, 3))
      d$y <- 2 * d$x + rnorm(n, 0, 0.5)
      list(f = bf(y ~ a * (exp(b * x) - 1), a + b ~ 1, nl = TRUE), d = d,
           fam = gaussian(), pr = NULL)
    },
    fn_hillstep = {
      d <- data.frame(x = c(runif(n / 2, 0.5, 4.5), runif(n / 2, 5.5, 10)))
      d$y <- 10 * (d$x > 5) + rnorm(n, 0, 0.5)
      list(f = bf(y ~ Emax * x^h / (EC50^h + x^h), Emax + EC50 + h ~ 1,
                  nl = TRUE), d = d, fam = gaussian(), pr = NULL)
    },
    fn_mmlin = {
      d <- data.frame(x = runif(n, 0, 1))
      d$y <- 0.5 * d$x + rnorm(n, 0, 0.2)
      list(f = bf(y ~ Vm * x / (K + x), Vm + K ~ 1, nl = TRUE), d = d,
           fam = gaussian(), pr = NULL)
    },
    fn_poismmlin = {
      d <- data.frame(x = runif(n, 0.1, 5))
      d$y <- rpois(n, 0.8 * d$x)
      list(f = bf(y ~ log(Vm * x / (K + x)), Vm + K ~ 1, nl = TRUE),
           d = d, fam = poisson(), pr = NULL)
    })
}
# REVIEW_PLACE = prior: a located normal(truth, 2 |truth| + 1) prior on
# each nlpar the model sets none on; start: start at the truth
place <- Sys.getenv("REVIEW_PLACE", "none")
truth <- list(logis = c(Asym = 10, xmid = 5, scal = 1),
              poismm = c(Vm = 8, K = 1.5), bernlogis = c(m = 4, s = 1.2),
              hill = c(Emax = 10, EC50 = 3, h = 2),
              fn_explin = c(a = 2, b = 0.5),
              fn_hillstep = c(Emax = 10, EC50 = 5, h = 2),
              fn_mmlin = c(Vm = 1, K = 1), fn_poismmlin = c(Vm = 1, K = 1))
cat("place:", place, "\n")
cat("ROWHEAD model n seed secs ll code wclass mclass est\n")
for (m in models) for (n in ns_) for (seed in seeds) {
  g <- gen(m, n, seed)
  st <- NULL
  if (place == "prior" && !is.null(truth[[m]])) {
    tr <- truth[[m]]
    g$pr <- do.call(c, lapply(names(tr), function(p) {
      set_prior(sprintf("normal(%g, %g)", tr[[p]], 2 * abs(tr[[p]]) + 1),
                nlpar = p)
    }))
  } else if (place == "start" && !is.null(truth[[m]])) {
    st <- list(beta = unname(truth[[m]]))
  }
  t0 <- proc.time()[["elapsed"]]
  r <- tt(frm(g$f, family = g$fam, data = g$d, prior = g$pr, start = st))
  secs <- proc.time()[["elapsed"]] - t0
  mc <- attr(r, "mc")
  est <- if (iserr(r)) r$msg else
    paste(format(unname(r$estimates$beta), digits = 5), collapse = ",")
  cat(sprintf("ROW %s %d %d %.2f %.6f %s %s %s %s\n", m, n, seed, secs,
              llv(r), if (iserr(r)) "ERR" else r$opt$convergence, wsum(r),
              if (length(mc)) paste(unique(mc), collapse = "+") else "-",
              est))
  if (any(attr(r, "wc") == "frmtmb_nl_limit")) {
    cat("  LIMIT:", attr(r, "w")[attr(r, "wc") == "frmtmb_nl_limit"], "\n")
  }
  flush.console()
}
cat("DONE\n")
