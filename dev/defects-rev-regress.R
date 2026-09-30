# Reviewer of lane defects: base against lane on many model shapes.
# Each arm fits the same models and saves every output; a third call
# compares them with identical().
#   Rscript dev/defects-rev-regress.R lane|base   (writes the rds)
#   Rscript dev/defects-rev-regress.R compare
arm <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- "dev/defects-rev-log"
if (identical(arm, "compare")) {
  a <- readRDS(file.path(out_dir, "regress-base.rds"))
  b <- readRDS(file.path(out_dir, "regress-lane.rds"))
  cat("models:", length(a), "\n")
  nd <- 0; ni <- 0
  for (m in names(a)) {
    for (k in union(names(a[[m]]), names(b[[m]]))) {
      ni <- ni + 1
      x <- a[[m]][[k]]; y <- b[[m]][[k]]
      if (!identical(x, y)) {
        nd <- nd + 1
        ae <- tryCatch(isTRUE(all.equal(x, y, check.attributes = FALSE)),
                       error = function(e) FALSE)
        cat(sprintf("DIFF %-14s %-18s all.equal(no attr)=%s\n", m, k, ae))
        if (is.character(x) || is.character(y)) {
          dx <- setdiff(x, y); dy <- setdiff(y, x)
          if (length(dx)) cat("   base only:", head(dx, 6), sep = "\n     ")
          if (length(dy)) cat("   lane only:", head(dy, 6), sep = "\n     ")
          cat("\n")
        } else if (!ae) {
          cat("   base:", substr(paste(format(x), collapse = " "), 1, 300), "\n")
          cat("   lane:", substr(paste(format(y), collapse = " "), 1, 300), "\n")
        } else {
          cat("   attributes base:", paste(names(attributes(x)), collapse = ","),
              " lane:", paste(names(attributes(y)), collapse = ","), "\n")
          if (!identical(dimnames(x), dimnames(y))) {
            cat("   dimnames differ\n")
          }
        }
      }
    }
  }
  cat("items compared:", ni, " different:", nd, "\n")
  quit(save = "no")
}
libs <- c("C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-defects-lib", libs)
.libPaths(libs)
suppressMessages(library(frmtmb))
set.seed(4242)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), w = runif(n),
                g = factor(rep(1:12, each = 10)),
                h = factor(rep(1:8, 15)),
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                lg = sample(c(TRUE, FALSE), n, TRUE),
                tt = rep(1:10, 12))
re_g <- rnorm(12, 0, 0.7)[d$g]
d$y <- 1 + 0.5 * d$x - 0.3 * d$z + re_g + rnorm(n, 0, 0.8)
d$y2 <- 0.5 * d$y + rnorm(n)
d$cnt <- rpois(n, exp(0.5 + 0.3 * d$x))
d$zcnt <- ifelse(runif(n) < 0.3, 0L, d$cnt)
d$tr <- 10L
d$bin <- rbinom(n, d$tr, plogis(0.2 + 0.5 * d$x))
d$yb <- rbinom(n, 1, plogis(0.5 * d$x))
d$pos <- exp(0.3 + 0.4 * d$x + rnorm(n, 0, 0.5))
d$hpos <- ifelse(runif(n) < 0.25, 0, d$pos)
d$prop <- plogis(0.2 * d$x + rnorm(n, 0, 0.5))
d$ord <- cut(d$y, quantile(d$y, c(0, 0.25, 0.5, 0.75, 1)),
             include.lowest = TRUE, labels = FALSE)
d$ordf <- factor(d$ord, ordered = TRUE)
d$cat <- factor(sample(c("p", "q", "r"), n, TRUE))
d$catn <- sample(c(2L, 5L, 9L), n, TRUE)
d$mixy <- ifelse(runif(n) < 0.5, rnorm(n, -2, 0.7), rnorm(n, 2, 0.7))
d$sei <- runif(n, 0.2, 0.6)
d$sdx <- 0.2
d$xm <- d$x + rnorm(n, 0, 0.2)
d$xmi <- ifelse(runif(n) < 0.15, NA, d$x)
d$wt <- runif(n, 0.5, 2)
d$off <- log(runif(n, 1, 3))
d$mof <- factor(sample(1:4, n, TRUE), ordered = TRUE)
d2 <- d
contrasts(d2$f) <- contr.sum(3)
models <- list(
  gauss = quote(frm(y ~ x + z, data = d)),
  gauss_f = quote(frm(y ~ f + x, data = d)),
  gauss_lg = quote(frm(y ~ lg + x, data = d)),
  gauss_contr = quote(frm(y ~ f + x, data = d2)),
  re1 = quote(frm(y ~ x + (1 | g), data = d)),
  re2 = quote(frm(y ~ x + (1 + x | g) + (1 | h), data = d)),
  dist = quote(frm(bf(y ~ x, sigma ~ z), data = d)),
  pois = quote(frm(cnt ~ x, family = poisson(), data = d)),
  pois_re = quote(frm(cnt ~ x + (1 | g), family = poisson(), data = d)),
  binom = quote(frm(bin | trials(tr) ~ x, family = binomial(), data = d)),
  bern = quote(frm(yb ~ x, family = bernoulli(), data = d)),
  negbin = quote(frm(cnt ~ x, family = negbinomial(), data = d)),
  lognorm = quote(frm(pos ~ x, family = lognormal(), data = d)),
  student = quote(frm(y ~ x, family = student(), data = d)),
  beta = quote(frm(prop ~ x, family = Beta(), data = d)),
  cumul = quote(frm(ordf ~ x + z, family = cumulative(), data = d)),
  sratio_cs = quote(frm(ord ~ x + cs(z), family = sratio(), data = d)),
  acat_cs = quote(frm(ord ~ x + cs(z), family = acat(), data = d)),
  cratio = quote(frm(ord ~ x + (1 | g), family = cratio(), data = d)),
  categ = quote(frm(cat ~ x, family = categorical(), data = d)),
  mixture = quote(frm(mixy ~ 1, family = mixture(gaussian(), gaussian()),
                      data = d)),
  mv_rescor = quote(frm(bf(y ~ x) + bf(y2 ~ z) + set_rescor(TRUE), data = d)),
  mv_norescor = quote(frm(bf(y ~ x) + bf(cnt ~ z, family = poisson()) +
                            set_rescor(FALSE), data = d)),
  gp = quote(frm(y ~ gp(x), data = d)),
  smooth = quote(frm(y ~ s(x) + z, data = d)),
  smooth2 = quote(frm(y ~ s(x, z), data = d)),
  nonlin = quote(frm(bf(y ~ b0 + b1 * exp(b2 * x), b0 + b1 + b2 ~ 1,
                        nl = TRUE), data = d,
                     prior = c(set_prior("normal(1, 2)", nlpar = "b0"),
                               set_prior("normal(0, 2)", nlpar = "b1"),
                               set_prior("normal(0, 1)", nlpar = "b2")))),
  arma = quote(frm(y ~ x + arma(tt, g), data = d)),
  ar1 = quote(frm(y ~ x + ar(tt, g, p = 1), data = d)),
  me = quote(frm(y ~ me(xm, sdx) + z, data = d)),
  mi = quote(frm(bf(y ~ mi(xmi) + z) + bf(xmi | mi() ~ z) +
                   set_rescor(FALSE), data = d)),
  zip = quote(frm(bf(zcnt ~ x, zi ~ z), family = zero_inflated_poisson(),
                  data = d)),
  hurdle_p = quote(frm(zcnt ~ x, family = hurdle_poisson(), data = d)),
  hurdle_ln = quote(frm(hpos ~ x, family = hurdle_lognormal(), data = d)),
  se = quote(frm(y | se(sei) ~ x, data = d)),
  se_sigma = quote(frm(y | se(sei, sigma = TRUE) ~ x, data = d)),
  weights = quote(frm(y | weights(wt) ~ x, data = d)),
  mo = quote(frm(y ~ mo(mof) + x, data = d)),
  offset = quote(frm(cnt ~ x + offset(off), family = poisson(), data = d)),
  cs_re = quote(frm(ord ~ cs(x) + (1 | g), family = sratio(), data = d))
)
safe <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)),
           error = function(e) paste("ERROR:", conditionMessage(e)))
}
res <- list()
nd <- d[c(3, 17, 45, 88), ]
for (nm in names(models)) {
  cat(nm, "\n")
  r <- list()
  fit <- safe(eval(models[[nm]]))
  if (is.character(fit)) {
    res[[nm]] <- list(fit = fit)
    next
  }
  p <- fit$opt$par %||% fit$obj$par
  r$par <- p
  r$fn <- safe(fit$obj$fn(p))
  r$gr <- safe(fit$obj$gr(p))
  r$logLik <- safe(as.numeric(logLik(fit)))
  r$fixef <- safe(fixef(fit))
  r$fitted <- safe(fitted(fit))
  r$fitted_lin <- safe(fitted(fit, scale = "linear"))
  r$fitted_nd <- safe(fitted(fit, newdata = nd))
  r$fitted_slice <- safe(fitted(fit, newdata = fit$data[1:20, ]))
  r$predict <- safe({set.seed(1); predict(fit, ndraws = 50)})
  r$predict_nd <- safe({set.seed(2); predict(fit, newdata = nd, ndraws = 50)})
  r$resid <- safe(residuals(fit))
  r$resid_p <- safe(residuals(fit, type = "pearson"))
  r$ranef <- safe(ranef(fit))
  r$prior <- safe(as.data.frame(default_prior(fit$bform, data = d)))
  r$summary <- safe(capture.output(print(summary(fit))))
  r$print <- safe(capture.output(print(fit)))
  r$family <- safe(capture.output(print(family(fit))))
  r$vc <- safe(capture.output(print(VarCorr(fit))))
  r$nobs <- safe(nobs(fit))
  res[[nm]] <- r
}
saveRDS(res, file.path(out_dir, paste0("regress-", arm, ".rds")))
cat("done", length(res), "models\n")
