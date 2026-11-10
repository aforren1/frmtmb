# Reviewer probe: does draws_newdata_is_fitted() accept newdata that is
# not the fitted rows, or refuse the true fitted data?
.libPaths(c("C:/Users/adf44/source/r/wt-hollow-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb.sample"), "\n")
isf <- frmtmb.sample:::draws_newdata_is_fitted
mk <- function(fit) structure(list(stanfit = NULL, draws = matrix(0, 1, 1),
                                   fit = fit), class = "frmtmb_draws")
q <- function(expr) suppressWarnings(suppressMessages(expr))

set.seed(1)
n <- 40
dd <- data.frame(
  x = rnorm(n), z = rnorm(n), w = runif(n, 0.5, 2),
  g = factor(rep(letters[1:5], each = 8)), t = rep(1:8, 5),
  nn = rep(10L, n), s = rep(c(TRUE, TRUE, TRUE, FALSE), 10),
  x1 = runif(n), x2 = runif(n), sdx = runif(n, 0.1, 0.3),
  age = runif(n, 1, 3), cc = rep(c(0, 0, 0, 2), 10), u = runif(n))
dd$y <- 1 + 0.5 * dd$x + rnorm(n)
dd$y2 <- dd$y + 1
dd$yc <- dd$y
dd$k <- rbinom(n, 10, 0.4)
dd$ypos <- exp(dd$y / 3)
dd$xo <- factor(sample(1:4, n, TRUE), ordered = TRUE)

perturb <- function(d, v) {
  a <- d[[v]]
  if (identical(v, "t")) {
    # swap two time points within the first group
    a[c(1, 2)] <- a[c(2, 1)]
  } else if (is.factor(a)) {
    a[c(1, n)] <- a[c(n, 1)]
  } else if (is.logical(a)) {
    a[1] <- !a[1]
  } else if (is.integer(a)) {
    a[1] <- a[1] + 1L
  } else {
    a[2] <- a[2] + 0.25
  }
  d[[v]] <- a
  d
}

cases <- list(
  ar = function() frm(y ~ x + ar(time = t, gr = g, p = 1), data = dd),
  arma_noTime = function() frm(y ~ x + ar(gr = g, p = 1), data = dd),
  cosy = function() frm(y ~ x + cosy(time = t, gr = g), data = dd),
  weights = function() frm(y | weights(w) ~ x, data = dd),
  trials = function() frm(k | trials(nn) ~ x, family = binomial(), data = dd),
  cens = function() frm(y | cens(cc, y2) ~ x, data = dd),
  se = function() frm(y | se(sdx) ~ x, data = dd),
  subset = function() frm(bf(y | subset(s) ~ x), data = dd),
  gp = function() frm(y ~ gp(x1, x2), data = dd),
  me = function() frm(y ~ me(x, sdx), data = dd),
  nl = function() frm(bf(y ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
                     data = dd),
  nl_cov = function() frm(bf(y ~ a * exp(b * x), a ~ 1 + z, b ~ 1,
                             nl = TRUE), data = dd),
  dpar = function() frm(bf(y ~ x, sigma ~ z), data = dd),
  poly_off = function() frm(y ~ poly(x, 2) + offset(log(age)), data = dd),
  logy = function() frm(log(ypos) ~ x, data = dd),
  mo = function() frm(y ~ mo(xo), data = dd),
  mm = function() frm(y ~ x + (1 | mm(g, g2)),
                      data = transform(dd, g2 = rev(g))),
  grby = function() frm(y ~ x + (1 | gr(g, by = s)), data = dd),
  mv_subset = function() frm(bf(y | subset(s) ~ x) + bf(z ~ x) +
                               set_rescor(FALSE), data = dd)
)

res <- list()
for (nm in names(cases)) {
  fit <- tryCatch(q(cases[[nm]]()), error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-11s FIT ERROR: %s\n", nm, conditionMessage(fit)))
    next
  }
  ds <- mk(fit)
  dat <- if (nm == "mm") transform(dd, g2 = rev(g)) else dd
  mfn <- names(fit$data)
  rawn <- names(fit$frame$raw_vars)
  fast <- all(mfn %in% names(dat))
  cat(sprintf("%-11s mf=[%s] raw=[%s] data_vars=[%s] fast=%s\n", nm,
              paste(mfn, collapse = ","), paste(rawn, collapse = ","),
              paste(fit$frame$data_vars, collapse = ","), fast))
  cat(sprintf("   accepts data=%s fit$data=%s\n", isf(ds, dat),
              isf(ds, fit$data)))
  # a perturbation of every column of dat; flag acceptance of a change
  # to a variable the formula names
  for (v in names(dat)) {
    acc <- isf(ds, perturb(dat, v))
    if (acc) cat(sprintf("   ACCEPTS change to %s\n", v))
  }
}
