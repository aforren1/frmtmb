# Review of lane wcheck: assorted attacks.
#   Rscript dev/reviews/wcheck-review-misc.R <lib|base> <part>
# part jac:   natural-scale priors on sigma toward 0 (MAP carries the
#             log-Jacobian), verdict vs the profile slope
# part cens:  a bounded intercept prior that blocks the mean shift;
#             asym_laplace censored on both sides (errs open?)
# part quad:  what quadseam can still use as a quadrature runaway
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
          substr(m, 1, 150))
  } else {
    sprintf("fit logLik %.6g betad %s conv %s warn %d%s",
            as.numeric(logLik(r)),
            paste(format(r$estimates$betad, digits = 4), collapse = " "),
            r$opt$convergence, length(ws),
            if (length(ws)) paste0(" (", substr(ws[length(ws)], 1, 70), ")")
            else "")
  }
}
profile_disp <- function(ob, g, j = 1L, fix = NULL) {
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  isd <- which(names(p) == "betad")[j]
  if (!is.null(fix)) p[names(fix)] <- fix
  free <- setdiff(seq_along(p), c(isd, match(names(fix), names(p))))
  vapply(g, function(t) {
    p[isd] <- t
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
                g = factor(rep(1:16, each = 5)), t = rep(1:5, 16))
d$y <- 1 + 0.5 * d$x + rnorm(n)

if (part == "jac") {
  gs <- -c(2, 4, 6, 8, 10, 12)
  for (E in c(0.5, 0.9)) {
    d$w <- 1 + E / n
    fo <- bf(y | weights(w) ~ x + (1 | obs))
    for (ps in c("normal(0, 1)", "student_t(3, 0, 2.5)", "exponential(1)")) {
      pr <- set_prior(ps, class = "sigma")
      ob <- q(frm(fo, family = gaussian(), data = d, prior = pr,
                  dry_run = "objective"))
      show(sprintf("gaussian (1 | obs), E %g, %s on sigma | check: %s", E, ps,
                   outcome(frm(fo, family = gaussian(), data = d,
                               prior = pr))),
           -gs, profile_disp(ob, gs))
    }
  }
}

if (part == "cens") {
  q8 <- quantile(d$y, 0.8)
  d$cc <- ifelse(d$y > q8, "right", "none")
  d$yc <- pmin(d$y, q8)
  nu <- sum(d$cc == "none")
  E <- 2
  d$w <- ifelse(d$cc == "none", 1 + E / nu, 1)
  cat(sprintf("bound %.3f, %d right censored, E %g\n", q8, sum(d$cc != "none"),
              E))
  fo <- bf(yc | weights(w) + cens(cc) ~ 1 + (1 | obs))
  for (pr in list(NULL,
                  set_prior("uniform(-5, 0.5)", class = "Intercept",
                            lb = -5, ub = 0.5))) {
    lab <- if (is.null(pr)) "no prior" else "uniform(-5, 0.5) on Intercept"
    gs <- -c(2, 4, 6, 8, 10, 12, 14)
    ob <- tryCatch(q(frm(fo, family = gaussian(), data = d, prior = pr,
                         dry_run = "objective")), error = function(e) e)
    if (inherits(ob, "error")) {
      cat(lab, ": objective ERROR", conditionMessage(ob), "\n")
      next
    }
    show(sprintf("gaussian right-censored, y ~ 1, %s | check: %s", lab,
                 outcome(frm(fo, family = gaussian(), data = d, prior = pr))),
         -gs, profile_disp(ob, gs))
  }
  # asym_laplace with 36 right, 36 left censored and 8 uncensored at w 1:
  # the lane's cost is 36 halves, excess 8 - 18 < 0, let through
  o <- order(d$y)
  d$ca <- "none"
  d$ca[o[1:36]] <- "left"
  d$ca[o[45:80]] <- "right"
  lo <- d$y[o[36]]; hi <- d$y[o[45]]
  d$yac <- pmin(pmax(d$y, lo), hi)
  fa <- bf(yac | cens(ca) ~ 1 + (1 | obs))
  gs <- -c(2, 4, 6, 8, 10, 12, 14)
  ob <- q(frm(fa, family = asym_laplace(), data = d, dry_run = "objective"))
  show(sprintf("asym_laplace 36 left + 36 right, 8 uncensored | check: %s",
               outcome(frm(fa, family = asym_laplace(), data = d))),
       -gs, profile_disp(ob, gs))
  ob2 <- if (!is.null(ob$obj)) ob$obj else ob
  iq <- which(names(ob2$par) == "betad")[2L]
  fix <- setNames(0, names(ob2$par)[iq])
  # hold the quantile at the median by index
  pp <- ob2$par
  prof_q <- function(g) {
    p <- pp; isd <- which(names(p) == "betad")
    free <- setdiff(seq_along(p), isd)
    p[isd[2L]] <- 0
    vapply(g, function(t) {
      p[isd[1L]] <- t
      fr <- function(r) {
        q2 <- p; q2[free] <- r
        v <- tryCatch(ob2$fn(q2), error = function(e) NA_real_)
        if (!is.finite(v)) 1e10 else v
      }
      op <- nlminb(p[free], fr)
      if (op$objective < 1e9) p[free] <<- op$par
      -op$objective
    }, 0)
  }
  show("  same, quantile held at 0.5", -gs, prof_q(gs))
  cat("  frm(quadrature = TRUE):",
      outcome(frm(fa, family = asym_laplace(), data = d, quadrature = TRUE)),
      "\n")
}

if (part == "quad") {
  d$w <- 1.5
  d$s <- 0.3
  d$ym <- d$y; d$ym[c(5, 17, 33)] <- NA
  d$yz <- ifelse(runif(n) < 0.2, 0, exp(d$y / 2))
  d$yzb <- plogis(d$y - 1); d$yzb[1:10] <- 0
  cands <- list(
    list("se(s, sigma = TRUE)",
         bf(y | se(s, sigma = TRUE) + weights(w) ~ x + (1 | obs)),
         gaussian()),
    list("ar(time = t, gr = g)",
         bf(y | weights(w) ~ x + (1 | obs) + ar(time = t, gr = g)),
         gaussian()),
    list("mi() with 3 missing",
         bf(ym | weights(w) + mi() ~ x + (1 | obs)), gaussian()),
    list("hurdle_lognormal", bf(yz | weights(w) ~ x + (1 | obs)),
         hurdle_lognormal()),
    list("zero_inflated_beta", bf(yzb | weights(w) ~ x + (1 | obs)),
         zero_inflated_beta()))
  for (cs in cands) {
    for (quad in c(FALSE, TRUE)) {
      cat(sprintf("%-26s quadrature %-5s: %s\n", cs[[1]], quad,
                  outcome(frm(cs[[2]], family = cs[[3]], data = d,
                              quadrature = quad))))
    }
  }
}
cat("DONE\n")
