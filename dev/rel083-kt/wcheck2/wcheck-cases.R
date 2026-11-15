# Lane wcheck: the weights check's verdict on every case of items 3447,
# 3450/3608, 3435 and their guard-absent cases, before and after.
#   Rscript dev/wcheck-cases.R <lib|base> <part: mv|cens|prior>
# Each ARM line is the fit's outcome: REFUSED (the weights check), ERROR
# (anything else), or ok with logLik, the dispersion coefficients
# (betad), the nlminb code and the warnings. The data of each part are
# those of the record's scripts, with their seeds:
#   mv:    dev/fuzzfix-rev-wu.R (seed 20261021, n = 80)
#   cens:  dev/rel073-wcens-probe.R (seed 20261021, n = 80)
#   prior: dev/fuzzfix-rev-r2-edges.R and dev/fuzzfix-rev-r1-priors.R
#          (seed 20261021, n = 80)
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
part <- a[2]
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "| part", part, "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
run <- function(name, expr) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    m <- gsub("[[:space:]]+", " ", conditionMessage(r))
    tag <- if (grepl("has no maximum|is improper", m)) "REFUSED" else "ERROR"
    cat(sprintf("ARM %-46s %s | %s\n", name, tag, substr(m, 1, 160)))
  } else {
    cat(sprintf("ARM %-46s ok logLik %.6g | betad %s | conv %s | warnings %d%s\n",
                name, as.numeric(logLik(r)),
                paste(sprintf("%.4g", r$estimates$betad), collapse = " "),
                r$opt$convergence %||% NA, length(ws),
                if (length(ws)) paste0(" (", substr(ws[1], 1, 60), ")") else ""))
  }
  invisible(r)
}
G <- gaussian()
P <- function(...) set_prior(...)

if (part == "mv") {
  set.seed(20261021)
  n <- 80
  d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n), xs = rnorm(n),
                  w = runif(n, 0.5, 2), g = factor(rep(1:20, each = 4)))
  eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
  d$y <- eta + rnorm(n, 0, exp(0.3 * d$xs))
  d$yg <- rgamma(n, shape = 3, rate = 3 / exp(eta))
  d$k <- rpois(n, exp(eta))
  d$y2 <- d$y + rnorm(n)
  # an overdispersed count: a row-level gaussian effect on the log mean
  d$k2 <- rpois(n, exp(eta + rnorm(n, 0, 0.6)))
  d$k3 <- rpois(n, exp(0.3 + 0.2 * d$x + rnorm(n, 0, 0.6)))
  d$w7 <- d$w * 0.7
  d$sub <- seq_len(n) <= 60
  cat(sprintf("sum(w) - n = %.4f; sum(w7) - n = %.4f; rows 1-60 %.4f\n",
              sum(d$w) - n, sum(d$w7) - n, sum(d$w[d$sub] - 1)))
  cat("--- runaways: the check must refuse\n")
  run("mv gaussian x2 (1 | obs), w",
      frm(mvbf(bf(y | weights(w) ~ x + (1 | obs)),
               bf(y2 | weights(w) ~ x + (1 | obs))), family = G, data = d))
  run("mv gaussian + poisson, w",
      frm(mvbf(bf(y | weights(w) ~ x + (1 | obs)) + G,
               bf(k | weights(w) ~ x + (1 | obs)) + poisson()), data = d))
  run("mv poisson + Gamma, w",
      frm(mvbf(bf(k | weights(w) ~ x + (1 | obs)) + poisson(),
               bf(yg | weights(w) ~ x + (1 | obs)) + Gamma("log")),
          data = d))
  run("mv gaussian x2 (1 | p | obs), w",
      frm(mvbf(bf(y | weights(w) ~ x + (1 | p | obs)),
               bf(y2 | weights(w) ~ x + (1 | p | obs))), family = G,
          data = d))
  run("mv gaussian x2, w on y2 only",
      frm(mvbf(bf(y ~ x + (1 | obs)),
               bf(y2 | weights(w) ~ x + (1 | obs))), family = G, data = d))
  run("mv gaussian, subset(sub) on y, w",
      frm(mvbf(bf(y | weights(w) + subset(sub) ~ x + (1 | obs)),
               bf(y2 | weights(w7) ~ x + (1 | obs))), family = G, data = d))
  run("mv gaussian x2, lognormal prior sigma resp y",
      frm(mvbf(bf(y | weights(w) ~ x + (1 | obs)),
               bf(y2 | weights(w) ~ x + (1 | obs))), family = G, data = d,
          prior = P("lognormal(0, 1)", class = "sigma", resp = "y")))
  run("mv gaussian x2 (1 | obs), w, quadrature",
      frm(mvbf(bf(y | weights(w) ~ x + (1 | obs)),
               bf(y2 | weights(w) ~ x + (1 | obs))), family = G, data = d,
          quadrature = TRUE))
  cat("--- a maximum exists: the check must not refuse\n")
  run("mv gaussian x2 (1 | obs), w7 (sum < n)",
      frm(mvbf(bf(y | weights(w7) ~ x + (1 | obs)),
               bf(y2 | weights(w7) ~ x + (1 | obs))), family = G, data = d))
  run("mv gaussian x2 (1 | g), w",
      frm(mvbf(bf(y | weights(w) ~ x + (1 | g)),
               bf(y2 | weights(w) ~ x + (1 | g))), family = G, data = d))
  run("mv gaussian x2, lognormal prior sigma both",
      frm(mvbf(bf(y | weights(w) ~ x + (1 | obs)),
               bf(y2 | weights(w) ~ x + (1 | obs))), family = G, data = d,
          prior = c(P("lognormal(0, 1)", class = "sigma", resp = "y"),
                    P("lognormal(0, 1)", class = "sigma", resp = "y2"))))
  run("poisson (1 | obs), overdispersed, unweighted",
      frm(k2 ~ x + (1 | obs), family = poisson(), data = d))
  run("poisson (1 | obs), overdispersed, w",
      frm(k2 | weights(w) ~ x + (1 | obs), family = poisson(), data = d))
  run("mv poisson x2 (1 | obs), overdispersed",
      frm(mvbf(bf(k2 ~ x + (1 | obs)), bf(k3 ~ x + (1 | obs))),
          family = poisson(), data = d))
  run("mv poisson x2 (1 | obs), overdispersed, w",
      frm(mvbf(bf(k2 | weights(w) ~ x + (1 | obs)),
               bf(k3 | weights(w) ~ x + (1 | obs))),
          family = poisson(), data = d))
  run("mv poisson x2 (1 | p | obs), overdispersed",
      frm(mvbf(bf(k2 ~ x + (1 | p | obs)), bf(k3 ~ x + (1 | p | obs))),
          family = poisson(), data = d))
  run("mv gaussian + poisson (1 | obs), unweighted",
      frm(mvbf(bf(y ~ x + (1 | g)) + G,
               bf(k2 ~ x + (1 | obs)) + poisson()), data = d))
  run("mv gaussian x2 (1 | obs), unweighted",
      frm(mvbf(bf(y ~ x + (1 | obs)), bf(y2 ~ x + (1 | obs))),
          family = G, data = d))
  run("mv gaussian x2 (1 | obs), w7, quadrature",
      frm(mvbf(bf(y | weights(w7) ~ x + (1 | obs)),
               bf(y2 | weights(w7) ~ x + (1 | obs))), family = G, data = d,
          quadrature = TRUE))
  m <- tryCatch(frm(mvbf(bf(y ~ x + (1 | obs)),
                         bf(y2 | weights(w) ~ x + (1 | obs))),
                    family = G, data = d),
                error = function(e) conditionMessage(e),
                warning = function(w) conditionMessage(w))
  cat("MESSAGE mv:", m, "\n")
}

if (part == "cens") {
  set.seed(20261021)
  n <- 80
  d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)),
                  g4 = factor(rep(1:20, each = 4)))
  d$y <- 1 + 0.5 * d$x + rnorm(n)
  q <- quantile(d$y, 0.8)
  d$cc <- ifelse(d$y > q, "right", "none")
  d$yc <- pmin(d$y, q)
  nu <- sum(d$cc == "none")
  cat("--- gaussian, 16 of 80 right censored, uncensored excess E\n")
  for (E in c(0, 0.5, 2, 4, 7.5, 8, 16)) {
    d$w <- ifelse(d$cc == "none", 1 + E / nu, 1)
    run(sprintf("gaussian right E %g Laplace", E),
        frm(bf(yc | weights(w) + cens(cc) ~ x + (1 | obs)),
            family = G, data = d))
  }
  d$w <- ifelse(d$cc == "none", 1 + 4 / nu, 1)
  run("gaussian right E 4 quadrature",
      frm(bf(yc | weights(w) + cens(cc) ~ x + (1 | obs)),
          family = G, data = d, quadrature = TRUE))
  run("gaussian right E 4 (1 | g4)",
      frm(bf(yc | weights(w) + cens(cc) ~ x + (1 | g4)),
          family = G, data = d))
  # both sides: the 8 lowest left censored too, so no one shift of the
  # mean puts every censored row's mass at 1; the allowance stays 1/2
  # per row on the smaller side
  ql <- quantile(d$y, 0.1)
  d$c2 <- ifelse(d$y > q, "right", ifelse(d$y < ql, "left", "none"))
  d$y2c <- pmin(pmax(d$y, ql), q)
  nu2 <- sum(d$c2 == "none")
  cat(sprintf("both sides: %d right, %d left, %d uncensored\n",
              sum(d$c2 == "right"), sum(d$c2 == "left"), nu2))
  for (E in c(2, 3.5, 4.5, 8)) {
    d$w2 <- ifelse(d$c2 == "none", 1 + E / nu2, 1)
    run(sprintf("gaussian both E %g Laplace", E),
        frm(bf(y2c | weights(w2) + cens(c2) ~ x + (1 | obs)),
            family = G, data = d))
  }
  cat("--- asym_laplace, unweighted, k of 80 right censored\n")
  d$ya <- 1 + 0.5 * d$x + rnorm(n)
  for (k in c(0, 8, 40, 72, 80)) {
    o <- order(d$ya, decreasing = TRUE)
    ca <- rep("none", n)
    if (k) ca[o[seq_len(k)]] <- "right"
    qa <- if (k) d$ya[o[k]] else Inf
    d$ca <- ca
    d$yac <- pmin(d$ya, qa)
    run(sprintf("asym_laplace %2d of 80 right censored", k),
        frm(yac | cens(ca) ~ x + (1 | obs), family = asym_laplace(),
            data = d))
  }
  run("asym_laplace 72 right censored (1 | g4)",
      frm(yac | cens(ca) ~ x + (1 | g4), family = asym_laplace(), data = d))
}

if (part == "prior") {
  set.seed(20261021)
  n <- 80
  d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                  fac = factor(rep(c("a", "b"), 40)),
                  w = runif(n, 0.5, 2), g = factor(rep(1:20, each = 4)))
  eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
  d$y <- eta + rnorm(n, 0, 0.8)
  d$yg <- rgamma(n, shape = 3, rate = 3 / exp(eta))
  d$yb <- plogis(eta); d$yb <- rbeta(n, d$yb * 5, (1 - d$yb) * 5)
  cat(sprintf("excess: all %.3f | fac a %.3f | fac b %.3f\n", sum(d$w) - n,
              sum(d$w[d$fac == "a"] - 1), sum(d$w[d$fac == "b"] - 1)))
  fo <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 0 + fac)
  cat("--- cell means: sigma ~ 0 + fac\n")
  run("cellmeans, no prior", frm(fo, family = G, data = d))
  run("cellmeans, normal(0,1) on b", frm(fo, family = G, data = d,
      prior = P("normal(0, 1)", class = "b", dpar = "sigma")))
  run("cellmeans, normal on faca only", frm(fo, family = G, data = d,
      prior = P("normal(0, 1)", class = "b", dpar = "sigma", coef = "faca")))
  run("cellmeans, normal on facb only", frm(fo, family = G, data = d,
      prior = P("normal(0, 1)", class = "b", dpar = "sigma", coef = "facb")))
  run("cellmeans, student_t on b", frm(fo, family = G, data = d,
      prior = P("student_t(3, 0, 2.5)", class = "b", dpar = "sigma")))
  run("cellmeans, normal faca + student_t facb", frm(fo, family = G, data = d,
      prior = c(P("normal(0, 1)", class = "b", dpar = "sigma", coef = "faca"),
                P("student_t(3, 0, 2.5)", class = "b", dpar = "sigma",
                  coef = "facb"))))
  cat("--- sigma ~ 0 + fac + x\n")
  fox <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 0 + fac + x)
  run("fac + x, normal(0,1) on b", frm(fox, family = G, data = d,
      prior = P("normal(0, 1)", class = "b", dpar = "sigma")))
  run("fac + x, normal on faca + facb", frm(fox, family = G, data = d,
      prior = c(P("normal(0, 1)", class = "b", dpar = "sigma", coef = "faca"),
                P("normal(0, 1)", class = "b", dpar = "sigma", coef = "facb"))))
  cat("--- sigma ~ 1 + fac\n")
  fo1 <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
  run("1 + fac, normal on Intercept only", frm(fo1, family = G, data = d,
      prior = P("normal(0, 1)", class = "Intercept", dpar = "sigma")))
  run("1 + fac, normal on Intercept and b", frm(fo1, family = G, data = d,
      prior = c(P("normal(0, 1)", class = "Intercept", dpar = "sigma"),
                P("normal(0, 1)", class = "b", dpar = "sigma"))))
  cat("--- Gamma shape, Beta phi (natural scale, toward infinity)\n")
  fg <- bf(yg | weights(w) ~ x + (1 | obs))
  GA <- Gamma("log")
  run("Gamma, student_t(3,0,2.5) shape", frm(fg, family = GA, data = d,
      prior = P("student_t(3, 0, 2.5)", class = "shape")))
  run("Gamma, student_t(30,0,2.5) shape", frm(fg, family = GA, data = d,
      prior = P("student_t(30, 0, 2.5)", class = "shape")))
  run("Gamma, inv_gamma(2,1) shape", frm(fg, family = GA, data = d,
      prior = P("inv_gamma(2, 1)", class = "shape")))
  run("Gamma, gamma(0.01,0.01) shape", frm(fg, family = GA, data = d,
      prior = P("gamma(0.01, 0.01)", class = "shape")))
  fb <- bf(yb | weights(w) ~ x + (1 | obs))
  run("Beta, student_t(3,0,2.5) phi", frm(fb, family = Beta(), data = d,
      prior = P("student_t(3, 0, 2.5)", class = "phi")))
  run("Beta, gamma(0.01,0.01) phi", frm(fb, family = Beta(), data = d,
      prior = P("gamma(0.01, 0.01)", class = "phi")))
  cat("--- Gamma shape ~ 0 + fac\n")
  fgc <- bf(yg | weights(w) ~ x + (1 | obs), shape ~ 0 + fac)
  run("Gamma cellmeans, normal(0,1) on b", frm(fgc, family = GA, data = d,
      prior = P("normal(0, 1)", class = "b", dpar = "shape")))
  run("Gamma cellmeans, normal on faca only", frm(fgc, family = GA, data = d,
      prior = P("normal(0, 1)", class = "b", dpar = "shape", coef = "faca")))
  # backlog 3486's arm: a positive-support prior on the log-scale cells
  run("Gamma cellmeans, gamma(2,1) on b", frm(fgc, family = GA, data = d,
      prior = P("gamma(2, 1)", class = "b", dpar = "shape")))
  run("Gamma cellmeans, gamma(2,1) on b, unweighted",
      frm(bf(yg ~ x + (1 | obs), shape ~ 0 + fac), family = GA, data = d,
          prior = P("gamma(2, 1)", class = "b", dpar = "shape")))
}
cat("DONE\n")
