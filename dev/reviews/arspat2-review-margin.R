# Review of lane arspat2, attack 1/2: the margin route on cases the lane
# did not run. The mi() fit's -obj$fn against the fit on the observed rows
# with fcor(M_oo), at three outer points (lane's comparison). Data seed 4
# unless said. Base when REVLIB="".
source("dev/reviews/arspat2-review-common.R")
n <- 60
M <- 0.5^abs(outer(1:n, 1:n, "-"))
mk <- function(re_sd = 1, seed = 4, M = get("M", globalenv()), df = 4,
               sig_sd = 0) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), h = factor(rep(1:6, each = 10)),
                  k = factor(rep(1:5, length.out = n)), z = runif(n))
  d$y <- 1 + 0.5 * d$x + rnorm(6, 0, re_sd)[d$h] + rnorm(5)[d$k] +
    exp(sig_sd * rnorm(6)[d$h] + 0.3 * d$x) *
    as.vector(t(chol(M)) %*% rnorm(n)) * sqrt(df / rchisq(1, df))
  d
}
run <- function(lab, f, d, miss, fam = student(), REML = FALSE,
                profile = FALSE, Mm = M, nus = c(3, 6, 21), sigf = NULL,
                prior = NULL) {
  d$y[miss] <- NA
  ok <- !is.na(d$y)
  fo <- f
  fo[[3]] <- do.call(substitute, list(fo[[3]], list(M = quote(Mo))))
  fm <- update(f, y | mi() ~ .)
  bfm <- if (is.null(sigf)) bf(fm) else bf(fm, sigf)
  bfo <- if (is.null(sigf)) bf(fo) else bf(fo, sigf)
  ctl <- frmtmb_control(profile = profile)
  t0 <- proc.time()[["elapsed"]]
  a <- tryCatch(quiet(frm(bfm, data = d, data2 = list(M = Mm), family = fam,
                          REML = REML, control = ctl, prior = prior)),
                error = function(e) e)
  if (inherits(a, "error")) {
    cat(sprintf("%-40s mi fit ERROR: %s\n", lab, conditionMessage(a)))
    return(invisible())
  }
  ta <- proc.time()[["elapsed"]] - t0
  b <- quiet(frm(bfo, data = d[ok, ], data2 = list(Mo = Mm[ok, ok]),
                 family = fam, REML = REML, control = ctl, prior = prior))
  lab_a <- ns$outer_par_labels(a)
  same <- identical(names(a$opt$par), names(b$opt$par))
  inu <- if (identical(fam$family, "student")) grep("^nu", lab_a) else integer(0)
  tm <- ns$autocor_seam_t_missing(a$frame)
  route <- if (!length(tm)) "none" else if (is.list(tm$y)) {
    if (isTRUE(tm$y$margin)) "margin" else "constant"
  } else "constant(base)"
  df <- if (same) cmp_pts(a, b, inu, nus) else NA
  cat(sprintf("%-40s route %-14s m %2d  diff %s  ll mi %.6f obs %.6f  exact %s  %.1fs\n",
              lab, route, sum(!ok),
              paste(sprintf("%+.2e", df), collapse = " "),
              as.numeric(logLik(a)), as.numeric(logLik(b)),
              ns$mi_laplace_exact(a$frame), ta))
  invisible(list(a = a, b = b))
}
pr <- set_prior("gamma(2, 0.1)", class = "nu")
d <- mk()
set.seed(40); half <- sort(sample(n, n / 2))
run("m = n/2, (1|h)", y ~ x + (1 | h) + fcor(M), d, half)
run("clustered: group 1 all missing + 5", y ~ x + (1 | h) + fcor(M), d,
    c(1:10, 21, 33, 45, 50, 58))
run("nu near 1 (1.02, 1.1, 1.5)", y ~ x + (1 | h) + fcor(M),
    mk(df = 1.2), c(3, 9, 17, 25, 31, 44, 52), nus = c(1.02, 1.1, 1.5))
run("two RE terms (1|h) + (1|k)", y ~ x + (1 | h) + (1 | k) + fcor(M), d,
    c(2, 8, 15, 22, 29, 36, 43, 50, 57))
run("(x|h)", y ~ x + (x | h) + fcor(M), d, c(2, 8, 15, 22, 29, 36, 43, 50))
Mn <- 0.995^abs(outer(1:n, 1:n, "-"))
cat("kappa(near-singular M):", format(kappa(Mn, exact = TRUE), digits = 3), "\n")
run("near-singular M 0.995^|i-j|, (1|h)", y ~ x + (1 | h) + fcor(M),
    mk(M = Mn), c(2, 8, 15, 22, 29, 36, 43, 50), Mm = Mn)
run("profile, no RE", y ~ x + fcor(M), d, c(2, 8, 15, 22, 29, 36), profile = TRUE)
run("student REML, sigma ~ x", y ~ x + fcor(M), mk(), c(2, 8, 15, 22, 29, 36),
    REML = TRUE, sigf = sigma ~ x)
# gaussian: integrated value in a predictor other than mu
run("gaussian REML, sigma ~ x", y ~ x + fcor(M), mk(), c(2, 8, 15, 22, 29, 36),
    fam = gaussian(), REML = TRUE, sigf = sigma ~ x)
run("gaussian profile, sigma ~ x", y ~ x + fcor(M), mk(), c(2, 8, 15, 22, 29, 36),
    fam = gaussian(), profile = TRUE, sigf = sigma ~ x)
run("gaussian sigma ~ (1|p|h), mu (1|p|h)", y ~ x + (1 | p | h) + fcor(M),
    mk(sig_sd = 0.6), c(2, 8, 15, 22, 29, 36), fam = gaussian(),
    sigf = sigma ~ 1 + (1 | p | h))
run("gaussian sigma ~ (1|k), mu (1|h)", y ~ x + (1 | h) + fcor(M),
    mk(sig_sd = 0.6), c(2, 8, 15, 22, 29, 36), fam = gaussian(),
    sigf = sigma ~ 1 + (1 | k))
run("gaussian sigma ~ 1 + (1|h) REML", y ~ x + fcor(M),
    mk(sig_sd = 0.6), c(2, 8, 15, 22, 29, 36), fam = gaussian(), REML = TRUE,
    sigf = sigma ~ 1 + (1 | h))
# smooth / gp in sigma: are they integrated, and does the route see them?
for (sf in list(sigma ~ s(z, k = 5), sigma ~ gp(z, k = 5))) {
  dd <- mk(sig_sd = 0.6)
  dd$y[c(2, 8, 15, 22, 29, 36)] <- NA
  a <- tryCatch(quiet(frm(bf(y | mi() ~ x + fcor(M), sf), data = dd,
                          data2 = list(M = M), family = gaussian())),
                error = function(e) e)
  if (inherits(a, "error")) {
    cat(deparse(sf), "ERROR", conditionMessage(a), "\n")
    next
  }
  tm <- ns$autocor_seam_t_missing(a$frame)
  cat(sprintf("gaussian %s: random %s, re_blocks %d, route %s, exact %s\n",
              deparse(sf), paste(unique(names(a$obj$env$par)[a$obj$env$random]),
                                 collapse = ","),
              length(a$frame$re_blocks),
              if (length(tm)) tm$y$margin else "none",
              ns$mi_laplace_exact(a$frame)))
}
