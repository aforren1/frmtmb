# Fix round 2, items 2 and 5. Usage: REV_ARM=... [REV_MUTLIB=...] Rscript ... <mode>
# units: Student at y x 1e-6, data seeds 1:4, against the fit in units (F6).
# family: a family's own parameter with a natural scale far below 1e-4 (item 2).
# floor: ordinal y ~ (5 z + m), seeds 1:8, m 1e4 to 3e5, the silent part (item 2).
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
arm <- Sys.getenv("REV_ARM", "after")
set_libs(arm)
suppressPackageStartupMessages(library(frmtmb))
cat("lib", find.package("frmtmb"), "\n")
ses <- function(f) tryCatch(sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed)), error = function(e) NA)
fitw <- function(...) {
  w <- character(0)
  f <- withCallingHandlers(tryCatch(frm(...), error = function(e) e),
                           warning = function(x) {
                             w <<- c(w, substr(conditionMessage(x), 1, 60))
                             invokeRestart("muffleWarning")
                           })
  list(f = f, w = w)
}
mode <- commandArgs(TRUE)[[1]]
if (mode == "units") {
  for (s in 1:4) {
    set.seed(s)
    g <- factor(rep(seq_len(15L), each = 10))
    z <- rnorm(150)
    y1 <- 0.3 + 0.5 * z + rnorm(15, 0, 0.7)[g] + 0.5 * rt(150, 5)
    a <- fitw(bf(y ~ z + (1 | g), nu = 4), family = student(), data = data.frame(y = y1 * 1e-6, z, g), importance = 300L)
    r <- fitw(bf(y ~ z + (1 | g), nu = 4), family = student(), data = data.frame(y = y1, z, g), importance = 300L)
    if (inherits(a$f, "error")) {
      cat(sprintf("student y x 1e-6 s %d STOP: %s\n", s, substr(conditionMessage(a$f), 1, 130)))
    } else {
      cat(sprintf("student y x 1e-6 s %d rounds %d slope SE %.4f theta SE %.4f warns %d\n", s, a$f$importance$rounds,
                  ses(a$f)[[2]] * 1e6 / ses(r$f)[[2]], ses(a$f)[[length(ses(a$f))]] / ses(r$f)[[length(ses(r$f))]], length(a$w)))
    }
  }
}
if (mode == "family") {
  # the same model at K = 1 and K = 1e6: kappa' = kappa K, so a step of 1e-4 on kappa
  # moves the predictor by 100 at K = 1e6 and by 1e-4 at K = 1
  fam_k <- function(K) {
    frmtmb_family(
      paste0("tpois", format(K)), dpars = "mu", links = list(mu = "log"),
      lpdf = function(y, dpars, aterms, extra) {
        e <- log(dpars[["mu"]]) + extra$kappa * K
        y * e - exp(e) - lgamma(y + 1)
      },
      extra_pars = function(y, aterms) list(kappa = 0))
  }
  # a scale: sigma = exp(kappa K), gaussian with a fixed mu link
  fam_s <- function(K) {
    frmtmb_family(
      paste0("tgau", format(K)), dpars = "mu", links = list(mu = "identity"),
      lpdf = function(y, dpars, aterms, extra) {
        ls <- extra$kappa * K
        -0.5 * ((y - dpars[["mu"]]) / exp(ls))^2 - ls - 0.5 * log(2 * pi)
      },
      extra_pars = function(y, aterms) list(kappa = 0))
  }
  for (s in 1:4) {
    set.seed(s)
    g <- factor(rep(seq_len(15L), each = 10))
    x <- rnorm(150)
    w <- rbinom(150, 1, 0.5)
    y <- rpois(150, exp(0.5 + 0.4 * x + rnorm(15, 0, 0.7)[g]))
    for (K in c(1e3, 1e5, 1e6, 1e8)) {
      r <- fitw(bf(y ~ 0 + x + (1 | g)), family = fam_k(1), data = data.frame(y, x, g), importance = 300L)
      a <- fitw(bf(y ~ 0 + x + (1 | g)), family = fam_k(K), data = data.frame(y, x, g), importance = 300L)
      if (inherits(a$f, "error") || inherits(r$f, "error")) {
        cat(sprintf("pois kappa K %g s %d STOP: %s\n", K, s, substr(conditionMessage(if (inherits(a$f, "error")) a$f else r$f), 1, 150)))
        next
      }
      sa <- ses(a$f); sr <- ses(r$f)
      nm <- names(a$f$opt$par)
      ik <- which(nm == "kappa")
      cat(sprintf("pois kappa K %g s %d rounds %d capped %s conv %d grad %.3g dll %+.5f mcse %.4f SE slope %.4f kappa %.4f theta %.4f warns %d %s\n", K, s,
                  a$f$importance$rounds, a$f$importance$capped, a$f$opt$convergence, a$f$importance$grad,
                  as.numeric(logLik(a$f)) - as.numeric(logLik(r$f)), a$f$importance$mcse, sa[[1]] / sr[[1]],
                  sa[[ik]] * K / sr[[ik]], sa[[which(nm == "theta")]] / sr[[which(nm == "theta")]], length(a$w), paste(a$w, collapse = ";")))
    }
    yy <- 0.3 + 0.5 * x + rnorm(15, 0, 0.7)[g] + rnorm(150, 0, 0.5)
    for (K in c(1e3, 1e6)) {
      r <- fitw(bf(y ~ x + (1 | g)), family = fam_s(1), data = data.frame(y = yy, x, g), importance = 64L)
      a <- fitw(bf(y ~ x + (1 | g)), family = fam_s(K), data = data.frame(y = yy, x, g), importance = 64L)
      if (inherits(a$f, "error") || inherits(r$f, "error")) {
        cat(sprintf("gauss kappa K %g s %d STOP: %s\n", K, s, substr(conditionMessage(if (inherits(a$f, "error")) a$f else r$f), 1, 150)))
        next
      }
      sa <- ses(a$f); sr <- ses(r$f)
      ik <- which(names(a$f$opt$par) == "kappa")
      cat(sprintf("gauss scale kappa K %g s %d rounds %d capped %s conv %d dll %+.2g SE slope %.4f kappa %.4f warns %d %s\n", K, s, a$f$importance$rounds,
                  a$f$importance$capped, a$f$opt$convergence, as.numeric(logLik(a$f)) - as.numeric(logLik(r$f)), sa[[2]] / sr[[2]], sa[[ik]] * K / sr[[ik]], length(a$w), paste(a$w, collapse = ";")))
    }
  }
}
if (mode == "floor") {
  for (m in c(1e4, 3e4, 1e5, 3e5)) for (s in 1:8) {
    set.seed(s)
    g <- factor(rep(seq_len(15L), each = 10))
    z <- rnorm(150)
    lat <- 0.7 * z + rnorm(15, 0, 0.7)[g] + rlogis(150)
    y <- ordered(as.integer(cut(lat, c(-Inf, -1, 0.5, Inf))))
    d <- data.frame(y, w = 5 * z + m, wc = 5 * z, g)
    a <- fitw(bf(y ~ w + (1 | g)), family = cumulative(), data = d, importance = 300L)
    r <- fitw(bf(y ~ wc + (1 | g)), family = cumulative(), data = d, importance = 300L)
    if (inherits(a$f, "error")) {
      cat(sprintf("m %-6g s %d STOP %s\n", m, s, substr(conditionMessage(a$f), 1, 120)))
      next
    }
    sa <- ses(a$f); sr <- ses(r$f)
    cat(sprintf("m %-6g s %d rounds %d conv %d log sd %.3f (c %.3f, z %+.3f) slope z %+.3f slopeSE %.4f ll-c %+.4f mcse %.4f warns %d %s\n", m, s,
                a$f$importance$rounds, a$f$opt$convergence, a$f$opt$par[[2]], r$f$opt$par[[2]],
                (a$f$opt$par[[2]] - r$f$opt$par[[2]]) / sr[[2]], (a$f$opt$par[[1]] - r$f$opt$par[[1]]) / sr[[1]], sa[[1]] / sr[[1]],
                as.numeric(logLik(a$f)) - as.numeric(logLik(r$f)), r$f$importance$mcse, length(a$w), paste(a$w, collapse = ";")))
  }
}
