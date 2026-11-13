# Short check of lane wcheck's fix round: attacks on the fixes.
#   Rscript dev/reviews/wcheck-short-attacks.R <lib|base> <part>
# part couple: other ways a prior reaches the moving dispersion direction
# part names:  response names whose brms spelling collides
# part shift:  priors on what the mean shift moves
# part ident:  a natural-scale prior on sigma with the identity link
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
part <- a[2]
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "| part", part, "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
outcome <- function(expr) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")},
    message = function(m) invokeRestart("muffleMessage")),
    error = function(e) e)
  if (inherits(r, "error")) {
    m <- gsub("\\s+", " ", conditionMessage(r))
    paste(if (grepl("has no maximum|is improper", m)) "REFUSED" else "ERROR",
          substr(m, 1, 170))
  } else {
    sprintf("fit logLik %.6g betad %s conv %s warn %d", as.numeric(logLik(r)),
            paste(format(r$estimates$betad, digits = 4), collapse = " "),
            r$opt$convergence, length(ws))
  }
}
# the objective along one parameter (by index among the "betad" or
# "beta" entries), the rest re-optimized
profile <- function(fo, family, prior, data, g, which = "betad", j = 1L) {
  ob <- q(frm(fo, family = family, data = data, prior = prior,
              dry_run = "objective"))
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  k <- which(names(p) == which)[j]
  free <- setdiff(seq_along(p), k)
  vapply(g, function(t) {
    p[k] <- t
    fr <- function(r) {
      q2 <- p; q2[free] <- r
      v <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
      if (!is.finite(v)) 1e10 else v
    }
    op <- tryCatch(nlminb(p[free], fr), error = function(e) NULL)
    if (is.null(op)) return(NA_real_)
    if (op$objective < 1e9) p[free] <<- op$par
    -op$objective
  }, 0)
}
show <- function(label, g, v) {
  cat(sprintf("%s\n  grid %s\n  -fn  %s\n  slope %s\n", label,
              paste(g, collapse = " "),
              paste(format(v, digits = 7), collapse = " "),
              paste(format(diff(v) / diff(g), digits = 3), collapse = " ")))
}
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                fac = factor(rep(c("a", "b"), 40)))
d$y <- 0.5 + 0.4 * d$x + rnorm(n, 0, ifelse(d$fac == "a", 0.6, 1.2))
G <- gaussian()
PS <- function(...) set_prior(..., dpar = "sigma")

if (part == "couple") {
  d$w <- 1 + 8 / n
  fx <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + x)
  # a density on the intercept at zero (class b, coef Intercept) under a
  # centered design: the constant direction moves b0, which it bounds
  pr <- PS("normal(0, 1)", class = "b", coef = "Intercept")
  cat("sigma ~ 1 + x, normal(0, 1) class b coef Intercept, E 8:",
      outcome(frm(fx, family = G, data = d, prior = pr)), "\n")
  gs <- c(0, -2, -4, -6, -8, -10, -12)
  v <- tryCatch(profile(fx, G, pr, d, gs), error = function(e) {
    cat("  profile ERROR", conditionMessage(e), "\n"); NULL})
  if (!is.null(v)) show("  profile in sigma's b0 (the intercept at zero)",
                        gs, v)
  # the same density as class Intercept, for comparison
  pr2 <- PS("normal(0, 1)", class = "Intercept")
  cat("sigma ~ 1 + x, normal(0, 1) class Intercept, E 8:",
      outcome(frm(fx, family = G, data = d, prior = pr2)), "\n")
  # class-wide b on sigma ~ 1 + fac, cells weighted 0.8 / 1.05
  d$w <- ifelse(d$fac == "a", 0.8, 1.05)
  f1 <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
  for (lab in c("normal(0, 1) class b", "student_t(3, 0, 2.5) class b",
                "normal(0, 1) Intercept + student_t(3, 0, 2.5) b",
                "student_t(3, 0, 2.5) Intercept + normal(0, 1) b")) {
    pr <- switch(lab,
      "normal(0, 1) class b" = PS("normal(0, 1)", class = "b"),
      "student_t(3, 0, 2.5) class b" = PS("student_t(3, 0, 2.5)", class = "b"),
      "normal(0, 1) Intercept + student_t(3, 0, 2.5) b" =
        c(PS("normal(0, 1)", class = "Intercept"),
          PS("student_t(3, 0, 2.5)", class = "b")),
      "student_t(3, 0, 2.5) Intercept + normal(0, 1) b" =
        c(PS("student_t(3, 0, 2.5)", class = "Intercept"),
          PS("normal(0, 1)", class = "b")))
    cat(sprintf("sigma ~ 1 + fac, %s: %s\n", lab,
                outcome(frm(f1, family = G, data = d, prior = pr))))
  }
  # sigma ~ 0 + fac with class b entries
  f0 <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 0 + fac)
  cat("sigma ~ 0 + fac, normal(0, 1) class b:",
      outcome(frm(f0, family = G, data = d,
                  prior = PS("normal(0, 1)", class = "b"))), "\n")
  cat("sigma ~ 0 + fac, normal on faca only:",
      outcome(frm(f0, family = G, data = d,
                  prior = PS("normal(0, 1)", class = "b", coef = "faca"))),
      "\n")
}

if (part == "names") {
  d$w <- runif(n, 0.5, 2)
  d$y_a <- d$y; d$ya <- d$y + rnorm(n); d$y.a <- d$y + rnorm(n)
  P <- function(...) set_prior(...)
  cat("mvbf(y_a, ya), normal sigma resp ya only (bounds ya, y_a free):",
      outcome(frm(mvbf(bf(y_a | weights(w) ~ x + (1 | obs)),
                       bf(ya | weights(w) ~ x + (1 | obs))),
                  family = G, data = d,
                  prior = P("lognormal(0, 1)", class = "sigma", resp = "ya")))
      , "\n")
  cat("mvbf(y_a, y.a), both lognormal resp raw:",
      outcome(frm(mvbf(bf(y_a | weights(w) ~ x + (1 | obs)),
                       bf(y.a | weights(w) ~ x + (1 | obs))),
                  family = G, data = d,
                  prior = c(P("lognormal(0, 1)", class = "sigma", resp = "y_a"),
                            P("lognormal(0, 1)", class = "sigma",
                              resp = "y.a")))), "\n")
}

if (part == "shift") {
  q8 <- quantile(d$y, 0.8)
  d$cc <- ifelse(d$y > q8, "right", "none")
  d$yc <- pmin(d$y, q8)
  nu <- sum(d$cc == "none")
  d$w <- ifelse(d$cc == "none", 1 + 2 / nu, 1)
  fo1 <- bf(yc | weights(w) + cens(cc) ~ 1 + (1 | obs))
  fox <- bf(yc | weights(w) + cens(cc) ~ x + (1 | obs))
  P <- function(...) set_prior(...)
  arms <- list(
    list("~ 1, lb/ub on Intercept without a density", fo1,
         P("", class = "Intercept", lb = -5, ub = 0.5)),
    list("~ 1, ub = 0.5 only, normal density", fo1,
         P("normal(0, 1)", class = "Intercept", ub = 0.5)),
    list("~ 1, lb = -5 only (up free)", fo1,
         P("normal(0, 1)", class = "Intercept", lb = -5)),
    list("~ x, uniform(0, 1) on the slope, Intercept free", fox,
         P("uniform(0, 1)", class = "b", lb = 0, ub = 1)),
    list("~ x, normal(0, 1) on the Intercept", fox,
         P("normal(0, 1)", class = "Intercept")))
  for (ar in arms) {
    cat(sprintf("%-50s %s\n", ar[[1]],
                outcome(frm(ar[[2]], family = G, data = d, prior = ar[[3]]))))
  }
  # the lane's construction on the slope arm: the mean held above the
  # bound with the slope at 0.5 (inside its prior), log sigma down
  ob <- q(frm(fox, family = G, data = d, prior = arms[[4]][[3]],
              dry_run = "objective"))
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  nm <- names(p)
  ib <- which(nm == "beta"); isd <- which(nm == "betad")
  free <- setdiff(seq_along(p), c(ib, isd))
  cat("  par:", nm, "\n")
  gs <- c(-2, -4, -6, -8, -10, -12)
  for (M in c(q8 + 3)) {
    p[ib] <- c(M, 0.5)
    v <- vapply(gs, function(t) {
      p[isd] <- t
      fr <- function(r) {q2 <- p; q2[free] <- r
        vv <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
        if (!is.finite(vv)) 1e10 else vv}
      op <- nlminb(p[free], fr)
      if (op$objective < 1e9) p[free] <<- op$par
      -op$objective
    }, 0)
    show(sprintf("  slope arm, Intercept %.2f, slope 0.5, path in log sigma",
                 M), -gs, v)
  }
}

if (part == "ident") {
  d$w <- 1 + 0.5 / n
  fo <- bf(y | weights(w) ~ x + (1 | obs))
  fam <- gaussian(link_sigma = "identity")
  pr <- set_prior("normal(0, 1)", class = "sigma")
  cat("identity link sigma, E 0.5, normal(0, 1) on sigma:",
      outcome(frm(fo, family = fam, data = d, prior = pr)), "\n")
  gs <- c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
  v <- tryCatch(profile(fo, fam, pr, d, gs), error = function(e) {
    cat("  profile ERROR", conditionMessage(e), "\n"); NULL})
  if (!is.null(v)) show("  profile in sigma (identity), per unit -log sigma",
                        -log(gs), v)
}
cat("DONE\n")
