# Reviewer of lane optima, claim 1(b)-(d): mo() beyond the study's one
# design. Each case is fitted twice in one process (determinism) on one
# arm; the two arms run as two processes on the same data.
#   Rscript dev/optima-rev-mo-cases.R base|lane [case regex]
arm <- commandArgs(TRUE)[1]
pick <- commandArgs(TRUE)[2]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb", find.package("frmtmb"), "\n")
ms <- tryCatch(utils::getFromNamespace("mo_simplex", "frmtmb"),
               error = function(e) function(z) {
                 x <- exp(c(0, z))
                 x / sum(x)
               })

mk <- function(seed, n = 300, K = 4) {
  set.seed(seed)
  d <- data.frame(x1 = sample(0:(K - 1), n, TRUE),
                  x2 = sample(0:3, n, TRUE),
                  z = rnorm(n), f = factor(sample(c("a", "b", "c"), n, TRUE)),
                  g = factor(sample(1:20, n, TRUE)))
  d
}
cases <- list()
add <- function(name, fun) cases[[name]] <<- fun
# two mo() terms, one with a zero middle step (a face)
add("two_mo", function(s) {
  d <- mk(s)
  d$y <- c(0, 1, 1, 2)[d$x1 + 1] + c(0, 0.5, 1.5, 1.6)[d$x2 + 1] +
    0.3 * d$z + rnorm(nrow(d))
  list(f = bf(y ~ mo(x1) + mo(x2) + z), fam = gaussian(), d = d)
})
# mo() times a factor
add("mo_x_factor", function(s) {
  d <- mk(s)
  d$y <- c(0, 1, 1, 2)[d$x1 + 1] * (1 + (d$f == "b")) + rnorm(nrow(d))
  list(f = bf(y ~ mo(x1) * f), fam = gaussian(), d = d)
})
# mo() in sigma
add("mo_in_sigma", function(s) {
  d <- mk(s)
  d$y <- 0.5 * d$z + rnorm(nrow(d), sd = exp(c(0, 0.4, 0.4, 0.8)[d$x1 + 1]))
  list(f = bf(y ~ z, sigma ~ mo(x1)), fam = gaussian(), d = d)
})
# mo() in an nlpar
add("mo_in_nlpar", function(s) {
  d <- mk(s)
  d$y <- 2 * exp(-0.1 * c(0, 1, 1, 3)[d$x1 + 1]) * (1 + 0.5 * d$z) +
    rnorm(nrow(d), sd = 0.3)
  list(f = bf(y ~ a * (1 + bz * z), a ~ mo(x1), bz ~ 1, nl = TRUE),
       fam = gaussian(), d = d)
})
# many levels: D = 9, steps with zeros
add("many_levels", function(s) {
  d <- mk(s, n = 500, K = 10)
  st <- c(0, 0.5, 0.5, 0.5, 1.2, 1.2, 2, 2, 2, 3)
  d$y <- st[d$x1 + 1] + rnorm(nrow(d))
  list(f = bf(y ~ mo(x1)), fam = gaussian(), d = d)
})
# no effect at all (b near 0), the weakly identified simplex
add("no_effect", function(s) {
  d <- mk(s)
  d$y <- 0.3 * d$z + rnorm(nrow(d))
  list(f = bf(y ~ mo(x1) + z), fam = gaussian(), d = d)
})
# random intercept
add("mo_ranef", function(s) {
  d <- mk(s)
  u <- rnorm(20, sd = 0.7)
  d$y <- c(0, 1, 1, 2)[d$x1 + 1] + u[d$g] + rnorm(nrow(d))
  list(f = bf(y ~ mo(x1) * z + (1 | g)), fam = gaussian(), d = d)
})
# mo() in a group-level term
add("mo_in_ranef", function(s) {
  d <- mk(s)
  u <- rnorm(20, sd = 0.3)
  d$y <- (1 + u[d$g]) * c(0, 1, 1, 2)[d$x1 + 1] + rnorm(nrow(d))
  list(f = bf(y ~ mo(x1) + (mo(x1) | g)), fam = gaussian(), d = d)
})
# me() beside mo()
add("mo_me", function(s) {
  d <- mk(s)
  d$xt <- rnorm(nrow(d))
  d$sdx <- 0.3
  d$xo <- d$xt + rnorm(nrow(d), sd = 0.3)
  d$y <- c(0, 1, 1, 2)[d$x1 + 1] + 0.5 * d$xt + rnorm(nrow(d))
  list(f = bf(y ~ mo(x1) + me(xo, sdx)), fam = gaussian(), d = d)
})
# mi() beside mo()
add("mo_mi", function(s) {
  d <- mk(s)
  d$w <- rnorm(nrow(d))
  d$y <- c(0, 1, 1, 2)[d$x1 + 1] + 0.5 * d$w + rnorm(nrow(d))
  d$w[sample(nrow(d), 30)] <- NA
  list(f = bf(y ~ mo(x1) + mi(w)) + bf(w | mi() ~ 1), fam = gaussian(),
       d = d)
})
# Poisson, binomial, cumulative
add("poisson", function(s) {
  d <- mk(s)
  d$y <- rpois(nrow(d), exp(0.2 + c(0, 0.6, 0.6, 1)[d$x1 + 1]))
  list(f = bf(y ~ mo(x1) * z), fam = poisson(), d = d)
})
add("bernoulli", function(s) {
  d <- mk(s)
  d$y <- rbinom(nrow(d), 1, plogis(-0.5 + c(0, 1, 1, 1.5)[d$x1 + 1]))
  list(f = bf(y ~ mo(x1) * z), fam = bernoulli(), d = d)
})
add("cumulative", function(s) {
  d <- mk(s)
  lat <- c(0, 1, 1, 2)[d$x1 + 1] + 0.4 * d$z + rlogis(nrow(d))
  d$y <- as.integer(cut(lat, c(-Inf, 0, 1, 2, Inf)))
  list(f = bf(y ~ mo(x1) * z), fam = cumulative(), d = d)
})
# negative effect with a face
add("negative", function(s) {
  d <- mk(s)
  d$y <- -c(0, 1, 1, 2)[d$x1 + 1] + rnorm(nrow(d))
  list(f = bf(y ~ mo(x1) * z), fam = gaussian(), d = d)
})
# a large model: 5000 rows, two mo() terms, an interaction, ranef
add("large", function(s) {
  d <- mk(s, n = 5000)
  d$g <- factor(sample(1:100, 5000, TRUE))
  u <- rnorm(100, sd = 0.5)
  d$y <- c(0, 0.2, 0.2, 0.5)[d$x1 + 1] * (1 + 0.1 * d$z) +
    c(0, 0, 0.3, 0.3)[d$x2 + 1] + u[d$g] + rnorm(5000)
  list(f = bf(y ~ mo(x1) * z + mo(x2) + (1 | g)), fam = gaussian(), d = d)
})

run1 <- function(cs, extra = list()) {
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers(
    do.call(frm, c(list(cs$f, family = cs$fam, data = cs$d), extra)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  list(fit = fit, w = w, t = proc.time()[["elapsed"]] - t0)
}
desc <- function(nm, s, r, r2, tag = "") {
  fit <- r$fit
  if (inherits(fit, "error")) {
    cat(sprintf("CASE %s%s seed %d ERROR %s\n", nm, tag, s,
                substr(gsub("\n", " ", conditionMessage(fit)), 1, 200)))
    return(invisible())
  }
  ll <- as.numeric(logLik(fit))
  zn <- grep("^zeta", names(fit$estimates), value = TRUE)
  mw <- vapply(zn, function(z) min(ms(fit$estimates[[z]])), 0)
  se <- tryCatch(sqrt(diag(vcov(fit, full = TRUE))), error = function(e) NA)
  det <- if (inherits(r2$fit, "error")) "ERR" else
    paste(identical(ll, as.numeric(logLik(r2$fit))),
          identical(unname(fit$opt$par), unname(r2$fit$opt$par)))
  cat(sprintf(paste0("CASE %s%s seed %d ll %.10f code %d evals %s ",
                     "search %s minw %s nonfinite_se %d warn %d t %.2f ",
                     "det %s\n"),
              nm, tag, s, ll, fit$opt$convergence,
              format(fit$opt$evals %||% NA),
              if (is.null(fit$opt$mo_search)) "-" else
                paste(format(fit$opt$mo_search, digits = 3), collapse = "/"),
              paste(format(mw, digits = 2), collapse = ","),
              sum(!is.finite(se)), length(r$w), r$t, det))
  for (m in unique(r$w)) cat("   W:", substr(gsub("\n", " ", m), 1, 160), "\n")
}
`%||%` <- function(a, b) if (is.null(a)) b else a
for (nm in names(cases)) {
  if (!is.na(pick) && !grepl(pick, nm)) next
  seeds <- if (nm == "large") 1:2 else 1:4
  for (s in seeds) {
    cs <- cases[[nm]](s)
    r <- run1(cs)
    r2 <- run1(cs)
    desc(nm, s, r, r2)
    if (nm %in% c("two_mo", "negative")) {
      ra <- run1(cs, list(control = frmtmb_control(profile = TRUE)))
      desc(nm, s, ra, ra, "+profile")
      ra <- run1(cs, list(REML = TRUE))
      desc(nm, s, ra, ra, "+REML")
      ra <- run1(cs, list(control = frmtmb_control(optimizer = "optim")))
      desc(nm, s, ra, ra, "+optim")
    }
  }
}
