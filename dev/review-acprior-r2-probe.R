# Reviewer, fix round 2: where the innovation reading of a sigma prior
# applies, and where it must not. Resolved entries (innov or not) and
# refusals; frm()'s MAP penalty against brms's density by hand.
# Usage: Rscript dev/review-acprior-r2-probe.R before|after
# Log: dev/review-acprior-log/r2-probe-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(6)
ng <- 12; k <- 6
d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                x = rnorm(ng * k))
e <- numeric(nrow(d))
for (i in seq_len(nrow(d))) {
  e[i] <- if (d$t[i] == 1) rnorm(1) else 0.6 * e[i - 1] + rnorm(1, 0, 0.8)
}
d$y1 <- 1 + 0.5 * d$x + e
d$y2 <- 2 - 0.3 * d$x + rnorm(nrow(d))
d$y <- d$y1
show <- function(lab, fo, pr, family = gaussian()) {
  cat("\n---", lab, "\n")
  tryCatch({
    fit <- frm(fo, data = d, prior = pr, family = family,
               dry_run = "objective")
    ri <- frmtmb:::resolve_prior_input(fit, pr)
    for (en in ri$entries) {
      cat("  entry", en$comp, paste(en$idx, collapse = ","), en$dist$kind,
          "scale", en$scale %||% "-", "innov",
          if (is.null(en$innov)) "no" else paste0("yes(thetaac ",
            paste(en$innov$idx, collapse = ","), ")"), "\n")
    }
    if (length(ri$lower)) cat("  lower:", paste(names(ri$lower),
                                                signif(ri$lower, 6)), "\n")
  }, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
fmv <- bf(y1 ~ x + ar(t, g, cov = TRUE)) + bf(y2 ~ x) + set_rescor(FALSE)
sp <- function(...) set_prior("student_t(3, 0, 2.5)", class = "sigma", ...)
show("mv: sigma resp y1 (ar response)", fmv, sp(resp = "y1"))
show("mv: sigma resp y2 (no autocor)", fmv, sp(resp = "y2"))
show("mv: sigma, no resp", fmv, sp())
show("mv: sigma resp y2 lb = 0.2 (must be taken)", fmv,
     set_prior("normal(1, 1)", class = "sigma", resp = "y2", lb = 0.2))
show("mv: sigma resp y1 lb = 0.2 (refused)", fmv,
     set_prior("normal(1, 1)", class = "sigma", resp = "y1", lb = 0.2))
show("ar cov: sigma lb = 0 (taken)", y ~ x + ar(t, g, cov = TRUE),
     set_prior("normal(1, 1)", class = "sigma", lb = 0))
show("ar cov: sigma ub = 3 (refused)", y ~ x + ar(t, g, cov = TRUE),
     set_prior("normal(1, 1)", class = "sigma", ub = 3))
show("ar cov = FALSE: sigma lb = 0.2 (taken, no innov)", y ~ x + ar(t, g),
     set_prior("normal(1, 1)", class = "sigma", lb = 0.2))
show("ar p = 2 cov: sigma lb = 0.2 (taken, no innov)",
     y ~ x + ar(t, g, p = 2, cov = TRUE),
     set_prior("normal(1, 1)", class = "sigma", lb = 0.2))
show("cosy: sigma lb = 0.2 (taken, no innov)", y ~ x + cosy(t, g),
     set_prior("normal(1, 1)", class = "sigma", lb = 0.2))
show("no autocor: sigma lb = 0.2", y ~ x,
     set_prior("normal(1, 1)", class = "sigma", lb = 0.2))
show("ar cov, identity-link sigma, lb = 0", y ~ x + ar(t, g, cov = TRUE),
     set_prior("normal(1, 1)", class = "sigma", lb = 0),
     family = brmsfamily("gaussian", link_sigma = "identity"))
show("ar cov, identity-link sigma, normal(1,1)", y ~ x + ar(t, g, cov = TRUE),
     set_prior("normal(1, 1)", class = "sigma"),
     family = brmsfamily("gaussian", link_sigma = "identity"))
show("ar cov + sigma ~ x: class sigma", bf(y ~ x + ar(t, g, cov = TRUE),
                                           sigma ~ x), sp())
show("ar cov + sigma ~ x: Intercept dpar sigma",
     bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x),
     set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma"))
show("ar cov + sigma ~ x: b dpar sigma",
     bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x),
     set_prior("normal(0, 1)", class = "b", dpar = "sigma"))
show("ar cov, student: nu prior (no innov)", y ~ x + ar(t, g, cov = TRUE),
     set_prior("gamma(2, 0.1)", class = "nu"), family = student())
show("ar cov, student: sigma prior (innov)", y ~ x + ar(t, g, cov = TRUE),
     sp(), family = student())

# frm()'s MAP penalty: difference of differences between two points
# against brms's density, half-t on sigma_i = sigma_m c, plus log c and
# the log link's Jacobian (log sigma_m)
cat("\n=== MAP penalty against brms's density by hand\n")
for (fo in list(y ~ x + ar(t, g, cov = TRUE), y ~ x + ma(t, g, cov = TRUE),
                y ~ x + arma(t, g, cov = TRUE))) {
  f1 <- frm(fo, data = d, prior = sp())
  f0 <- frm(fo, data = d)
  p1 <- f1$obj$env$last.par.best
  p2 <- p1
  p2[names(p2) == "thetaac"] <- p2[names(p2) == "thetaac"] * 0.5 + 0.1
  ls_i <- which(names(p1) == "betad")
  stopifnot(length(ls_i) == 1L)
  p2[ls_i] <- p2[ls_i] + 0.3
  pen <- function(p) f1$obj$fn(p) - f0$obj$fn(p)
  hand <- function(p) {
    th <- p[names(p) == "thetaac"]
    ac <- f1$frame$autocor[[1]]
    nat <- unlist(frmtmb:::autocor_natural(th, ac))
    ar <- if (ac$struct %in% c("ar", "arma")) nat[1] else 0
    ma <- if (ac$struct == "ma") nat[1] else if (ac$struct == "arma") nat[2] else 0
    c0 <- sqrt((1 - ar^2) / (1 + ma^2 + 2 * ar * ma))
    sm <- exp(p[ls_i])
    -(log(stats::dt(sm * c0 / 2.5, 3) / 2.5) + log(c0) + log(sm))
  }
  dd <- (pen(p2) - pen(p1)) - (hand(p2) - hand(p1))
  cat(sprintf("%-40s penalty change %.6f, hand %.6f, gap %.2e\n",
              deparse(fo), pen(p2) - pen(p1), hand(p2) - hand(p1), dd))
}
