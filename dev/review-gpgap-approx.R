# Review of lane gpgap, item 7: approximate gp() on UNSORTED positions, brms's
# own program with no rewriting, through brms_lp_check() (which asserts).
# Plain gp(x, k = ) has only the sorted row in test-brms-likelihood.R; this is
# the row that file lacks. Seed 8101. Also a short-lengthscale design for the
# claim that gp(x, k = 30) reaches brms's numbers, against the closed form.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-approx.R approx.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
H <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-stan.R")) {
  sys.source(file.path("tests/testthat", h), envir = H)
}
set.seed(8101)
n <- 80
d <- data.frame(x = sample(seq(0, 10, by = 0.25), n, TRUE),
                z = runif(n, 0, 5),
                f = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sin(d$x) + 0.5 * cos(d$z) + rnorm(n, 0, 0.3)
d$sx <- sort(d$x)
stopifnot(is.unsorted(d$x))
rows <- list(
  k10 = list(y ~ gp(x, k = 10)),
  k10_by = list(y ~ gp(x, by = f, k = 8)),
  k5_2d = list(y ~ gp(x, z, k = 5)),
  k5_2d_aniso = list(y ~ gp(x, z, k = 5, iso = FALSE)),
  k10_grF = list(y ~ gp(x, k = 10, gr = FALSE))
)
for (nm in names(rows)) {
  fo <- rows[[nm]][[1]]
  r <- tryCatch({
    fit <- suppressWarnings(frm(do.call(frmtmb::bf, list(fo)), family = gaussian(),
                                data = d))
    out <- with_reporter("silent", H$brms_lp_check(brms::bf(fo), gaussian(), d, fit,
                                                   joint = TRUE))
    sprintf("const %.3e grad %.2e tol %.1e conv %d", out$measured_const, out$max_grad,
            1e-6 * abs(out$ours), fit$opt$convergence)
  }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 150)))
  cat(sprintf("APPROX unsorted %-14s %s\n", nm, r))
}

# Is k = 30 enough? A field whose length scale is short against the range.
closed <- function(x, y, s0, nug) {
  u <- sort(unique(x)); D2 <- outer(u, u, "-")^2; J <- match(x, u)
  nll <- function(p) {
    s2 <- exp(2 * p[2])
    K <- s2 * exp(-D2 / (2 * exp(2 * p[3]))); diag(K) <- diag(K) + nug(s2)
    V <- K[J, J]; diag(V) <- diag(V) + exp(2 * p[4])
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(1e10)
    r <- backsolve(R, y - p[1], transpose = TRUE)
    0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * length(y) * log(2 * pi)
  }
  best <- NULL
  for (st in list(c(0, 0, -1, log(s0)), c(0, 0.7, -1.5, log(s0)),
                  c(0, 0.3, -0.5, log(s0)))) {
    o <- nlminb(st, nll, control = list(iter.max = 800, eval.max = 1200))
    if (is.null(best) || o$objective < best$objective) best <- o
  }
  best
}
set.seed(8102)
for (w in c(2, 5, 8, 12)) {
  n2 <- 200
  x <- round(runif(n2, 0, 6), 3)
  e <- rnorm(n2)
  for (s in c(1e-2, 1e-3)) {
    dd <- data.frame(x = x, y = sin(w * x) + s * e)
    b <- closed(x, dd$y, s, function(s2) 1e-12)
    a <- closed(x, dd$y, s, function(s2) 1e-6 * s2)
    cat(sprintf("KCHK w %g s %.0e | closed brms: sigma %.4e sdgp %.3f rho %.3f ll %.2f | frm nugget: sigma %.4e\n",
                w, s, exp(b$par[4]), exp(b$par[2]), exp(b$par[3]), -b$objective,
                exp(a$par[4])))
    for (k in c(15, 30, 60, 100)) {
      fk <- tryCatch(suppressWarnings(frm(frmtmb::bf(y ~ gp(x, k = k)), data = dd)),
                     error = function(e) NULL)
      if (is.null(fk)) {
        cat("  k", k, "ERROR\n")
        next
      }
      hk <- Filter(function(z) identical(z$covstruct, "hsgp"), fk$frame$re_blocks)[[1]]
      cat(sprintf("  k %3d: sigma %.4e (/brms %.2f) sdgp %.3f ll %.2f conv %d\n", k,
                  as.numeric(sigma(fk)), as.numeric(sigma(fk)) / exp(b$par[4]),
                  exp(fk$estimates$theta[hk$theta_idx][[1]]), -fk$opt$objective,
                  fk$opt$convergence))
    }
  }
}

# And on the trend designs, where sigma is ordinary and the field large.
set.seed(5501)
x <- round(runif(60, 0, 6), 3); e <- rnorm(60)
for (lab in c("lin", "qd", "expo")) {
  y <- switch(lab, lin = 5 * x, qd = 1.5 * x^2, expo = exp(1.2 * x)) + 0.3 * e
  dd <- data.frame(x = x, y = y)
  for (k in c(30, 100)) {
    fk <- tryCatch(suppressWarnings(frm(frmtmb::bf(y ~ gp(x, k = k)), data = dd)),
                   error = function(e) NULL)
    if (!is.null(fk)) {
      hk <- Filter(function(z) identical(z$covstruct, "hsgp"), fk$frame$re_blocks)[[1]]
      cat(sprintf("TRENDK %s k %d: sigma %.4f sdgp %.4g ll %.3f conv %d\n", lab, k,
                  as.numeric(sigma(fk)), exp(fk$estimates$theta[hk$theta_idx][[1]]),
                  -fk$opt$objective, fk$opt$convergence))
    }
  }
  # the honest spelling of a trend: the slope as a fixed effect
  fx <- suppressWarnings(frm(frmtmb::bf(y ~ x + gp(x)), data = dd))
  bk <- Filter(function(z) identical(z$covstruct, "gp"), fx$frame$re_blocks)[[1]]
  cat(sprintf("TRENDX %s y ~ x + gp(x): sigma %.4f sdgp %.4g ll %.3f\n", lab,
              as.numeric(sigma(fx)), exp(fx$estimates$theta[bk$theta_idx][[1]]),
              -fx$opt$objective))
}
