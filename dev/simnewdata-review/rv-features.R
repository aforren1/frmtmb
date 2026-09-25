# simulate(newdata = ) across model features. For each: the draws at
# newdata rows that are NOT the fitted rows, against the expected
# response fitted(newdata = ) (max |z| of the row means), and, where the
# construction allows, the identity with simulate() at the fitted rows.
#   Rscript dev/simnewdata-review/rv-features.R > .../log/features.txt
source("dev/simnewdata-review/rv-prelude.R")
set.seed(11)
n <- 300
d <- data.frame(x = rnorm(n), g = factor(rep(1:12, length.out = n)),
                t = runif(n, 0.5, 3), nt = sample(3:12, n, TRUE),
                w = runif(n, 0.5, 2))
eta <- 0.3 + 0.5 * d$x + rnorm(12, 0, 0.4)[d$g]
d$yp <- rpois(n, d$t * exp(eta))
d$yb <- rbinom(n, d$nt, plogis(eta))
d$yg <- eta + rnorm(n, 0, 0.5)
d$yz <- ifelse(runif(n) < plogis(-1 + d$x), 0, rpois(n, exp(eta)))
d$yh <- ifelse(runif(n) < 0.3, 0, 1 + rpois(n, exp(eta)))
d$yc <- factor(sample(c("a", "b", "c"), n, TRUE))
d$ycen <- pmin(d$yg, 1); d$cen <- as.integer(d$yg > 1)
d$lb <- -1
d$ytr <- ifelse(d$yg < -1, -1 + abs(rnorm(n, 0, .1)), d$yg)
d$se <- runif(n, 0.2, 1)
d$ymeta <- eta + rnorm(n, 0, d$se)
d$xm <- d$x; d$xm[sample(n, 30)] <- NA
d$sdx <- 0.2
nd <- data.frame(x = c(-2, 0, 1.5, 2.5), g = factor(c(1, 2, 3, 3)),
                 t = c(1, 5, 0.2, 2), nt = c(1, 20, 5, 50),
                 w = 1, lb = -1, se = c(0.1, 0.5, 1, 2), cen = 0L,
                 sdx = 0.2, xm = c(-2, 0, 1.5, 2.5))
nsim <- 4000

row_check <- function(label, f, family = NULL, newdata = nd, id = TRUE,
                      data = d, ...) {
  fit <- tryCatch(suppressWarnings(
    if (is.null(family)) frm(f, data = data, ...) else
      frm(f, family = family, data = data, ...)),
    error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-28s FIT ERROR: %s\n", label,
                substr(conditionMessage(fit), 1, 140)))
    return(invisible())
  }
  idm <- if (id) {
    a <- tryCatch(simulate(fit, nsim = 3, seed = 4),
                  error = function(e) e)
    b <- tryCatch(simulate(fit, nsim = 3, seed = 4, newdata = data),
                  error = function(e) e)
    if (inherits(a, "error") || inherits(b, "error")) {
      paste("id ERROR:", substr(conditionMessage(
        if (inherits(b, "error")) b else a), 1, 80))
    } else if (identical(a, b)) "id identical" else {
      am <- as.matrix(as.data.frame(lapply(a, as.numeric)))
      bm <- as.matrix(as.data.frame(lapply(b, as.numeric)))
      sprintf("id max rel diff %.2g", max(abs(am - bm) / pmax(abs(am), 1)))
    }
  } else "id skipped"
  s <- tryCatch(simulate(fit, nsim = nsim, seed = 2, newdata = newdata),
                error = function(e) e)
  if (inherits(s, "error")) {
    cat(sprintf("%-28s %s | SIM ERROR: %s\n", label, idm,
                substr(conditionMessage(s), 1, 160)))
    return(invisible(fit))
  }
  fv <- tryCatch(fitted(fit, newdata = newdata), error = function(e) e)
  if (inherits(fv, "error")) {
    cat(sprintf("%-28s %s | draws ok, fitted(newdata) ERROR: %s\n", label,
                idm, substr(conditionMessage(fv), 1, 120)))
    return(invisible(fit))
  }
  if (length(dim(fv)) == 3L) {
    # category probabilities: compare draw frequencies per category
    lv <- dimnames(fv)[[3]]
    ex <- fv[, "Estimate", , drop = TRUE]
    m <- as.matrix(as.data.frame(lapply(s, function(v) {
      if (is.factor(v)) as.integer(v) else v
    })))
    if (ncol(m) == nsim) {
      fr <- t(apply(m, 1, function(r) tabulate(r, length(lv)) / nsim))
      z <- (fr - ex) / sqrt(pmax(ex * (1 - ex), 1e-12) / nsim)
    } else {
      # a matrix response: mean count per category
      z <- NA
    }
    cat(sprintf("%-28s %s | category freq max|z| %.2f\n", label, idm,
                max(abs(z))))
    return(invisible(fit))
  }
  m <- as.matrix(as.data.frame(lapply(s, as.numeric)))
  mu <- fv[, "Estimate"]
  z <- (rowMeans(m) - mu) / (apply(m, 1, sd) / sqrt(nsim))
  cat(sprintf("%-28s %s | row means %s | E[y] %s | max|z| %.2f\n", label,
              idm, paste(format(rowMeans(m), digits = 3), collapse = " "),
              paste(format(mu, digits = 3), collapse = " "), max(abs(z))))
  invisible(fit)
}

row_check("gaussian (1|g)", bf(yg ~ x + (1 | g)))
row_check("gaussian REML", bf(yg ~ x + (1 | g)), REML = TRUE)
row_check("gaussian profile", bf(yg ~ x + (1 | g)),
          control = frmtmb_control(profile = TRUE))
row_check("poisson offset(log(t))", bf(yp ~ x + offset(log(t)) + (1 | g)),
          poisson())
row_check("poisson rate(t)", bf(yp | rate(t) ~ x + (1 | g)), poisson())
row_check("binomial trials(nt)", bf(yb | trials(nt) ~ x + (1 | g)),
          binomial())
row_check("gaussian weights(w)", bf(yg | weights(w) ~ x + (1 | g)))
row_check("gaussian cens(cen) latent", bf(ycen | cens(cen) ~ x))
row_check("gaussian trunc(lb = lb)", bf(ytr | trunc(lb = lb) ~ x))
row_check("gaussian se(se)", bf(ymeta | se(se) ~ x))
row_check("gaussian se(se, sigma=T)", bf(ymeta | se(se, sigma = TRUE) ~ x))
row_check("zi poisson zi ~ x", bf(yz ~ x + (1 | g), zi ~ x),
          zero_inflated_poisson())
row_check("hurdle poisson hu ~ x", bf(yh ~ x, hu ~ x), hurdle_poisson())
row_check("mo(x)", bf(yg ~ mo(xo)),
          data = transform(d, xo = cut(x, 4, labels = FALSE)),
          newdata = data.frame(xo = c(1, 2, 3, 4)))
row_check("me(x, sdx)", bf(yg ~ me(x, sdx)))
row_check("mi(xm) predictor", bf(yg ~ mi(xm)) + bf(xm | mi() ~ 1),
          id = FALSE)
row_check("categorical", bf(yc ~ x), categorical())
row_check("sigma ~ x + (1|g)", bf(yg ~ x, sigma ~ x + (1 | g)))
row_check("student", bf(yg ~ x + (1 | g)), student())
row_check("lognormal", bf(exp(yg) ~ x + (1 | g)), lognormal())
row_check("mixture 2 gauss", bf(yg ~ x), mixture(gaussian(), gaussian()))

# multinomial: trials from newdata
set.seed(3)
dm <- data.frame(x = rnorm(200), N = sample(5:15, 200, TRUE))
dm$Y <- t(sapply(seq_len(200), function(i) {
  stats::rmultinom(1, dm$N[i], exp(c(0, 0.5 * dm$x[i], -0.3)))
}))
colnames(dm$Y) <- c("a", "b", "c")
fm <- tryCatch(frm(bf(Y | trials(N) ~ x), family = multinomial(), data = dm),
               error = function(e) e)
if (inherits(fm, "error")) cat("multinomial FIT ERROR", conditionMessage(fm), "\n") else {
  ndm <- data.frame(x = c(-1, 1), N = c(2L, 100L))
  sm <- tryCatch(simulate(fm, nsim = 500, seed = 1, newdata = ndm),
                 error = function(e) e)
  if (inherits(sm, "error")) cat("multinomial SIM ERROR", conditionMessage(sm), "\n") else {
    tots <- vapply(sm, function(v) rowSums(as.matrix(v)), numeric(2))
    cat("multinomial newdata trials 2, 100: row totals of draws range",
        paste(range(tots[1, ]), collapse = "-"), "and",
        paste(range(tots[2, ]), collapse = "-"), "\n")
    mc <- Reduce(`+`, lapply(sm, as.matrix)) / 500
    fvm <- fitted(fm, newdata = ndm)
    cat("  mean counts:\n"); print(round(mc, 3))
    cat("  fitted(newdata):\n"); print(round(fvm[, "Estimate", ], 3))
  }
}

# cs() ordinal with a group term at newdata
d$yo <- factor(cut(d$x + rlogis(n), c(-Inf, -0.5, 0.5, Inf),
                   labels = FALSE), ordered = TRUE)
row_check("sratio cs(x)", bf(yo ~ cs(x)), sratio())
row_check("cumulative + (1|g)", bf(yo ~ x + (1 | g)), cumulative())

# student ar at newdata, where the residual is a multivariate t
dd <- expand.grid(time = 1:5, g = factor(1:40))
dd$y <- rnorm(nrow(dd)) + rep(rnorm(40), each = 5)
row_check("student ar(time, g)", bf(y ~ 1 + ar(time, gr = g, cov = TRUE)),
          student(), data = dd,
          newdata = data.frame(time = c(1, 2, 3), g = factor("z")))
# ar(cov = FALSE), brms's default: the lagged-residual form
row_check("ar(cov = FALSE)", bf(y ~ 1 + ar(time, gr = g)), data = dd,
          newdata = data.frame(time = c(1, 2, 3), g = factor("z")))
