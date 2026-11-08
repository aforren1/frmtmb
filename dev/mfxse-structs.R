# predictions() standard error against fitted()'s Est.Error across group-level
# structures, families and readers, one arm. A ratio of 1 means
# marginaleffects reads the uncertainty fitted() reads.
#   MFX_ARM=before|after Rscript dev/mfxse-structs.R
# Data seed 9101.
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
set.seed(9101)
n <- 360
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(1:12, n, TRUE)),
                h = factor(sample(1:6, n, TRUE)),
                f = factor(sample(c("a", "b"), n, TRUE)),
                g2 = factor(sample(1:12, n, TRUE)))
ug <- rnorm(12, 0, 0.7)
uh <- rnorm(6, 0, 0.5)
eta <- 0.3 + 0.4 * d$x + ug[d$g] + uh[d$h]
d$y <- eta + rnorm(n, 0, 0.6)
d$cnt <- rpois(n, exp(eta))
d$ys <- eta + rnorm(n, 0, exp(-0.5 + 0.4 * ug[d$g]))
d$zc <- ifelse(runif(n) < plogis(-1 + ug[d$g]), 0L, rpois(n, exp(eta)))
d$ord <- cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE)
d$y2 <- 0.5 * d$x + ug[d$g] + rnorm(n, 0, 0.6)
d$ynl <- 2 * exp(0.3 * d$x + 0.2 * ug[d$g]) + rnorm(n, 0, 0.3)
d$w1 <- 0.5
d$w2 <- 0.5
nd <- d[c(3, 70, 150, 220, 300), ]
specs <- list(
  crossed = list(bf(y ~ x + (1 | g) + (1 | h))),
  nested = list(bf(y ~ x + (1 | h / g))),
  uncor = list(bf(y ~ x + (x || g))),
  by = list(bf(y ~ x + (1 | gr(g, by = f)))),
  mm = list(bf(y ~ x + (1 | mm(g, g2, weights = cbind(w1, w2))))),
  t2 = list(bf(y ~ t2(x, z, k = 4))),
  hsgp = list(bf(y ~ gp(x, k = 10) + (1 | g))),
  id_sigma = list(bf(ys ~ x + (1 | p | g), sigma ~ (1 | p | g))),
  zip = list(bf(zc ~ x + (1 | g), zi ~ (1 | g)), family = zero_inflated_poisson()),
  nl = list(bf(ynl ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
            start = NULL),
  mv = list(bf(y ~ x + (1 | g)) + bf(y2 ~ x + (1 | g)) + set_rescor(FALSE)),
  cumul = list(bf(ord ~ x + (1 | g)), family = cumulative())
)
for (k in names(specs)) {
  s <- specs[[k]]
  args <- list(s[[1L]], data = d)
  if (!is.null(s$family)) args$family <- s$family
  fit <- tryCatch(quiet(do.call(frm, args)), error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-9s FITERR %s\n", k, conditionMessage(fit)))
    next
  }
  ws <- character(0)
  resp <- if (k == "mv") "y" else NULL
  p <- tryCatch(withCallingHandlers(
    predictions(fit, newdata = nd, resp = resp),
    warning = function(w) {
      ws <<- c(ws, substr(conditionMessage(w), 1, 60))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  ft <- tryCatch(quiet(fitted(fit, newdata = nd, resp = resp)),
                 error = function(e) e)
  if (inherits(p, "error") || inherits(ft, "error")) {
    cat(sprintf("%-9s ERR %s\n", k, conditionMessage(
      if (inherits(p, "error")) p else ft)))
    next
  }
  fse <- if (length(dim(ft)) == 3L) as.numeric(ft[, "Est.Error", ]) else
    ft[, "Est.Error"]
  r <- p$std.error / fse
  cat(sprintf("%-9s ncoef %4d  ratio %s  warn: %s\n", k,
              length(get_coef(fit)),
              paste(sprintf("%.4f", r), collapse = " "),
              paste(unique(ws), collapse = " | ")))
}
