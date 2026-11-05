# F1's second half and F2, F6: every class of outer parameter with no column,
# what step it now takes and what it reads, and large and small response units.
# Log dev/impsamp-out/classes.txt
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
ses <- function(f) sqrt(diag(ns$sdr_of(f)$cov.fixed))
fit1 <- function(form, fam, d, N = 300L, ...) {
  w <- character(0)
  f <- withCallingHandlers(
    tryCatch(suppressMessages(frm(form, family = fam, data = d, importance = N, ...)),
             error = function(e) e),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  list(fit = f, warn = w)
}
show <- function(lab, a, b, mult = 1) {
  if (inherits(b$fit, "error")) {
    cat(sprintf("%-46s STOP: %s\n", lab, substr(conditionMessage(b$fit), 1, 110)))
    return(invisible())
  }
  cls <- names(b$fit$opt$par)
  st <- ns$imp_outer_steps(b$fit$frame, b$fit$opt$par)
  sa <- ses(a$fit)
  sb <- ses(b$fit)
  cat(sprintf("%-46s rounds %d capped %s warns %d dll %+.5f (mcse %.4f)\n", lab,
              b$fit$importance$rounds, b$fit$importance$capped, length(b$warn),
              as.numeric(logLik(b$fit)) - as.numeric(logLik(a$fit)),
              a$fit$importance$mcse))
  for (cp in unique(cls)) {
    i <- which(cls == cp)
    cat(sprintf("    %-8s par %s step %s SE ratio %s\n", cp,
                paste(format(b$fit$opt$par[i], digits = 3), collapse = " "),
                paste(format(st[i], digits = 2), collapse = " "),
                paste(format(sb[i] * (if (cp == "beta") mult else 1) / sa[i],
                             digits = 4), collapse = " ")))
  }
}
ordd <- function(s, m) {
  set.seed(s)
  g <- factor(rep(seq_len(15L), each = 10))
  z <- rnorm(150)
  lat <- 0.7 * z + rnorm(15, 0, 0.7)[g] + rlogis(150)
  data.frame(y = ordered(as.integer(cut(lat, c(-Inf, -1, 0.5, Inf)))), w = 5 * z + m,
             g = g, h = factor(rep(c("a", "b"), 75)))
}
# ordinal thresholds, over families
for (fam in list(cumulative(), sratio(), cratio(), acat())) {
  a <- fit1(bf(y ~ w + (1 | g)), fam, ordd(5L, 0))
  for (m in c(1e4, 3e4)) {
    show(sprintf("%s, w + %g, data seed 5", fam$family, m), a,
         fit1(bf(y ~ w + (1 | g)), fam, ordd(5L, m)))
  }
}
# thres(gr =) holds several first thresholds
a <- fit1(bf(y | thres(gr = h) ~ w + (1 | g)), cumulative(), ordd(5L, 0))
if (inherits(a$fit, "error")) {
  cat("thres(gr = ):", substr(conditionMessage(a$fit), 1, 120), "\n")
} else {
  show("cumulative thres(gr = h), w + 1e4", a,
       fit1(bf(y | thres(gr = h) ~ w + (1 | g)), cumulative(), ordd(5L, 1e4)))
}
# theta at a correlation's limit and at a collapse
set.seed(69)
g <- factor(rep(seq_len(10L), each = 20))
x <- rnorm(200)
d <- data.frame(y = rpois(200, exp(0.4 + 0.5 * x + rnorm(10, 0, 0.5)[g] + 0.9 * rnorm(10, 0, 0.1)[g] * x)),
                x = x, g = g)
a <- fit1(bf(y ~ x + (1 + x | g)), poisson(), d, N = 2000L)
show("poisson (1 + x | g), theta: 200 against 2000 draws", a,
     fit1(bf(y ~ x + (1 + x | g)), poisson(), d))
# a slope's sd on a column of small and of large scale: theta's own size
set.seed(2)
g <- factor(rep(seq_len(20L), each = 8))
z <- rnorm(160)
d <- data.frame(y = rpois(160, exp(0.4 + 0.3 * z + rnorm(20, 0, 0.5)[g] + rnorm(20, 0, 0.3)[g] * z)),
                z = z, g = g)
a <- fit1(bf(y ~ z + (1 + z | g)), poisson(), d)
for (k in c(1e-6, 1e6)) {
  d2 <- d
  d2$z <- z * k
  show(sprintf("poisson (1 + z | g), z x %g: theta near %s", k, format(-log(k) - 1, digits = 3)),
       a, fit1(bf(y ~ z + (1 + z | g)), poisson(), d2), k)
}
# mo(): a simplex's chart
set.seed(3)
g <- factor(rep(seq_len(20L), each = 8))
o <- sample(0:3, 160, TRUE)
d <- data.frame(y = rpois(160, exp(0.2 + c(0, 0.3, 0.8, 1)[o + 1L] + rnorm(20, 0, 0.5)[g])), o = o, g = g)
a <- fit1(bf(y ~ mo(o) + (1 | g)), poisson(), d, N = 2000L)
show("poisson mo(o): 200 against 2000 draws", a, fit1(bf(y ~ mo(o) + (1 | g)), poisson(), d, N = 200L))
# large and small response units: gaussian and student, with no response unit in
# the step
for (fam in c("gaussian", "student")) for (k in c(1e-6, 1e-3, 1e3, 1e6, 1e7)) {
  res <- character(0)
  for (s in 1:4) {
    set.seed(s)
    g <- factor(rep(seq_len(15L), each = 10))
    z <- rnorm(150)
    u <- rnorm(15, 0, 0.7)[g]
    y0 <- 0.3 + 0.5 * z + u + (if (fam == "student") 0.5 * rt(150, 5) else rnorm(150, 0, 0.5))
    form <- if (fam == "student") bf(y ~ z + (1 | g), nu = 4) else bf(y ~ z + (1 | g))
    fm <- if (fam == "student") student() else gaussian()
    a <- fit1(form, fm, data.frame(y = y0, z = z, g = g), N = 200L)
    b <- fit1(form, fm, data.frame(y = y0 * k, z = z, g = g), N = 200L)
    if (inherits(b$fit, "error")) {
      res <- c(res, paste("STOP", substr(conditionMessage(b$fit), 1, 40)))
    } else {
      sb <- ses(b$fit)
      sa <- ses(a$fit)
      res <- c(res, sprintf("dll %+.1e SE %.4f", as.numeric(logLik(b$fit)) +
                              150 * log(k) - as.numeric(logLik(a$fit)),
                            sb[[2L]] / k / sa[[2L]]))
    }
  }
  cat(sprintf("%-8s y x %-6g %s\n", fam, k, paste(res, collapse = " | ")))
}
