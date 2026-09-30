# Reviewer: models that use none of the lane's four features must give
# identical() results on rellib-r3 and on the lane's build.
# Usage: Rscript dev/formula2-rev-regress.R before|after
arm <- commandArgs(TRUE)[1]
base <- c("C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "after") c("C:/Users/adf44/source/r/wt-formula2-lib",
                                base) else base)
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

set.seed(1)
n <- 240
g <- factor(rep(sprintf("g%02d", 1:12), each = n / 12))
f <- factor(sample(c("a", "b", "c"), n, TRUE))
f2 <- factor(sample(c("p", "q"), n, TRUE))
x <- rnorm(n); z <- rnorm(n)
ge <- rnorm(12, 0, 0.7)[as.integer(g)]
y <- 1 + 0.5 * x - 0.3 * z + c(0, 0.4, -0.6)[as.integer(f)] + ge +
  rnorm(n, 0, exp(0.2 * (f2 == "q")))
cnt <- rpois(n, exp(0.3 + 0.4 * x + ge / 2))
ord <- cut(y + rnorm(n), c(-Inf, 0, 1, 2, Inf), labels = FALSE)
cat3 <- factor(sample(c("u", "v", "w"), n, TRUE, prob = c(.5, .3, .2)))
t <- rep(1:(n / 12), 12)
sx <- runif(n, 0.1, 0.3); xo <- x + rnorm(n, 0, sx)
xm <- x; xm[sample(n, 20)] <- NA
comp <- rbinom(n, 1, 0.4)
ymix <- ifelse(comp == 1, 3 + 0.5 * x, -1 + 0.5 * x) + rnorm(n, 0, 0.8)
ynl <- 2 * exp(0.3 * x) + rnorm(n, 0, 0.3)
gby <- factor(ifelse(as.integer(g) %% 2 == 0, "e", "o"))
d <- data.frame(y, x, z, f, f2, g, cnt, ord, cat3, t, sx, xo, xm, ymix, ynl,
                gby)
nd <- d[c(1:10, 100:110), ]
nd_new <- nd; nd_new$g <- factor("gNEW")

grab <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e),
                                                 class = "ERR")),
    warning = function(w0) {
      w <<- c(w, conditionMessage(w0)); invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  list(value = v, warnings = w)
}

cases <- list(
  gauss_ix   = list(f = y ~ x * f + z),
  zero_x     = list(f = y ~ 0 + x),
  zero_icpt  = list(f = y ~ 0 + Intercept + x),
  cell_f     = list(f = y ~ 0 + f),
  cell_fx    = list(f = y ~ 0 + f:x),
  cell_ff    = list(f = y ~ 0 + f + f2),
  re_cell    = list(f = y ~ x + (0 + f | g)),
  re_slope   = list(f = y ~ x + (x | g)),
  sigma_f    = list(f = bf(y ~ x, sigma ~ 0 + f2)),
  nl         = list(f = bf(ynl ~ a * exp(b * x), a ~ 1, b ~ 1, nl = TRUE),
                    prior = c(prior(normal(2, 1), nlpar = "a"),
                              prior(normal(0, 1), nlpar = "b"))),
  mixture    = list(f = bf(ymix ~ x) + mixture(gaussian(), gaussian())),
  mixture_s  = list(f = bf(ymix ~ x, sigma2 ~ z) +
                      mixture(gaussian(), gaussian())),
  mixture_fx = list(f = bf(ymix ~ x, theta1 = 0.4) +
                      mixture(gaussian(), gaussian())),
  mv_plusfam = list(f = mvbf(bf(y ~ x), bf(cnt ~ x) + poisson()) +
                      gaussian()),
  mv_famarg  = list(f = mvbf(bf(y ~ x), bf(cnt ~ x) + poisson()),
                    family = gaussian()),
  mv_lf      = list(f = bf(y ~ x) + bf(z ~ x) +
                      lf(sigma ~ f2, resp = "y") + set_rescor(FALSE)),
  cs         = list(f = ord ~ cs(x) + z, family = acat()),
  thres      = list(f = ord | thres(3) ~ x, family = cumulative()),
  me         = list(f = y ~ me(xo, sx) + z),
  mi         = list(f = bf(y | mi() ~ mi(xm) + z) + bf(xm | mi() ~ z) +
                      set_rescor(FALSE)),
  gr_by      = list(f = y ~ x + (1 | gr(g, by = gby))),
  ar_nocov   = list(f = y ~ x + ar(t, g)),
  smooth     = list(f = y ~ s(x) + z),
  categ      = list(f = cat3 ~ x, family = categorical()),
  pois_re    = list(f = cnt ~ x + (1 | g) + offset(z / 10),
                    family = poisson()),
  lf_center  = list(f = bf(y ~ x) + lf(sigma ~ x, center = FALSE)),
  bf_center  = list(f = bf(y ~ x + f, center = FALSE)),
  bf_nested  = list(f = bf(bf(y ~ x), sigma ~ z))
)

out <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fam <- cs$family
  r <- list()
  fitr <- grab(frm(cs$f, data = d, family = fam, prior = cs$prior))
  fit <- fitr$value
  r$fit_warn <- fitr$warnings
  if (inherits(fit, "ERR")) {
    r$fit_err <- unclass(fit)
  } else {
    r$est <- grab(fit$estimates)
    r$opt <- grab(fit$opt$par)
    r$logLik <- grab(logLik(fit))
    r$fixef <- grab(fixef(fit))
    r$vcov <- grab(vcov(fit))
    r$variables <- grab(variables(fit))
    r$spec_pars <- grab(summary(fit)$spec_pars)
    r$lp_names <- names(fit$frame$linpreds)
    r$lp_idx <- lapply(fit$frame$linpreds, `[[`, "idx")
    r$lp_X <- lapply(fit$frame$linpreds, function(l) l[["X"]])
    r$blocks <- grab(lapply(fit$frame$re_blocks %||% fit$frame$blocks,
                            function(b) b$components))
    r$pred <- grab(predict(fit, newdata = nd))
    r$fitted <- grab(fitted(fit, newdata = nd))
    r$pred_new <- grab(predict(fit, newdata = nd_new,
                               allow_new_levels = TRUE))
    r$sim <- grab(simulate(fit, nsim = 2, seed = 7))
    r$ranef <- grab(ranef(fit))
    r$prior_fit <- grab(prior_summary(fit))
    r$upd_f <- grab(logLik(update(fit, . ~ . + f2)))
    r$upd_d <- grab(logLik(update(fit, newdata = d[-(1:12), ])))
  }
  r$default_prior <- grab(default_prior(cs$f, data = d, family = fam))
  r$par_template <- grab(par_template(cs$f, data = d, family = fam))
  np <- r$par_template$value
  r$frm_sim <- if (inherits(np, "ERR")) NULL else
    grab(frm_simulate(cs$f, data = d, family = fam, newparams = np,
                      seed = 3))
  out[[nm]] <- r
  cat(nm, if (!is.null(r$fit_err)) paste("ERR", r$fit_err) else "ok", "\n")
}
saveRDS(out, sprintf(
  "C:/Users/adf44/source/r/frmtmb-wt-formula2/dev/formula2-rev-regress-%s.rds",
  arm))
cat("DONE\n")
