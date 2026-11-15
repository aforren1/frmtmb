# Review of lane parity2, items 5315, 5300, 1560, 4160, 2980, 4337, 5323:
# frmtmb side (arm base|lane). brms side: parity2-review-names-brms.R.
#   Rscript dev/reviews/parity2-review-names.R base|lane
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-parity2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, " lib:", find.package("frmtmb"),
    as.character(packageVersion("frmtmb")), find.package("frmtmb.sample"),
    "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-parity2/dev/reviews/parity2-review-data.R")
try_ <- function(label, expr) {
  w <- character(); m <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    },
    message = function(cnd) {
      m <<- c(m, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    })
  if (inherits(r, "error")) cat(label, ": ERROR ", conditionMessage(r), "\n",
                                sep = "")
  for (x in unique(w)) cat("   warning: ", substr(x, 1, 160), "\n", sep = "")
  for (x in unique(m)) cat("   message: ", substr(gsub("\n", " ", x), 1, 160),
                           "\n", sep = "")
  invisible(r)
}
vars <- function(label, f, fam, data = rd, ...) {
  fit <- try_(label, frm(f, family = fam, data = data, ...))
  if (!inherits(fit, "error")) {
    cat(label, ": ", paste(grep("^(b_|bs|sd|sigma|Intercept|cor|r_)",
                                variables(fit), value = TRUE),
                           collapse = " "), "\n", sep = "")
  }
  invisible(fit)
}
rd$yo <- cut(rd$y, quantile(rd$y, c(0, 0.3, 0.6, 1)), include.lowest = TRUE,
             labels = FALSE)

cat("\n#### 5315 intercept place\n")
vars("center_false", bf(y ~ x, center = FALSE), gaussian())
vars("both_0I", bf(y ~ 0 + Intercept + x, sigma ~ 0 + Intercept + z),
     gaussian())
vars("ord_cf", bf(yo ~ x, center = FALSE), cumulative())
vars("mv_0I", bf(y ~ 0 + Intercept + x) + bf(y2 ~ z) + set_rescor(FALSE),
     gaussian())
vars("re_0I", bf(y ~ 0 + Intercept + x + (1 | g), sigma ~ z), gaussian())
vars("nl_sig0I", bf(y ~ a * exp(b * x), a ~ 1, b ~ 1,
                    sigma ~ 0 + Intercept + z, nl = TRUE), gaussian())

cat("\n#### 1560 nlf() replacing a parameter used in another's body\n")
f1 <- try_("nlf over lf a", bf(y ~ a + b * x, a ~ 1 + w, b ~ 1, nl = TRUE) +
             nlf(a ~ exp(c0)) + lf(c0 ~ 1 + z))
if (!inherits(f1, "error")) {
  fit <- try_("fit", frm(f1, family = gaussian(), data = rd))
  if (!inherits(fit, "error")) cat("variables:", paste(grep("^b_", variables(fit),
                                                       value = TRUE),
                                                  collapse = " "), "\n")
}
f2 <- try_("lf over fixed", bf(y ~ x, sigma = 1) + lf(sigma ~ z))
if (!inherits(f2, "error")) {
  cat("pfix sigma:", format(f2$pfix$sigma), " pforms sigma:",
      deparse1(f2$pforms$sigma), "\n")
}
f3 <- try_("lf over mu", bf(y ~ x) + lf(mu ~ z))
if (!inherits(f3, "error")) print(f3)

cat("\n#### 4160 attr(k, 'min') survives update() and drop1()\n")
set.seed(9)
d3 <- data.frame(k = rep(1:4, 10), z = rnorm(40))
d3$y <- 0.5 * d3$k + rnorm(40)
attr(d3$k, "min") <- 0L
f4 <- try_("mo attr", frm(bf(y ~ mo(k) + z), family = gaussian(), data = d3))
if (!inherits(f4, "error")) {
  m <- f4$frame$linpreds[[1]]$mo[[1]]
  cat("fit: min", m$min, " D", m$D, "\n")
  u <- try_("update", update(f4, . ~ . - z))
  if (!inherits(u, "error")) {
    m2 <- u$frame$linpreds[[1]]$mo[[1]]
    cat("update(. ~ . - z): min", m2$min, " D", m2$D, "\n")
  }
  u2 <- try_("update newdata", update(f4, newdata = d3[1:30, ]))
  if (!inherits(u2, "error")) {
    m3 <- u2$frame$linpreds[[1]]$mo[[1]]
    cat("update(newdata = d3[1:30, ]): min", m3$min, " D", m3$D, "\n")
  }
  dd <- try_("drop1", drop1(f4))
  if (!inherits(dd, "error")) print(dd)
  cat("logLik fit", as.numeric(logLik(f4)), "\n")
  d3b <- d3[sample(40), ]
  f4b <- try_("rows reordered", frm(bf(y ~ mo(k) + z), family = gaussian(),
                                    data = d3b))
  if (!inherits(f4b, "error")) {
    cat("reordered data min", f4b$frame$linpreds[[1]]$mo[[1]]$min, "\n")
  }
  d3n <- d3
  d3n$y[3] <- NA
  f4n <- try_("NA row", frm(bf(y ~ mo(k) + z), family = gaussian(),
                            data = d3n))
  if (!inherits(f4n, "error")) {
    cat("NA row min", f4n$frame$linpreds[[1]]$mo[[1]]$min, "\n")
  }
  d3t <- tibble::as_tibble(d3)
  f4t <- try_("tibble", frm(bf(y ~ mo(k) + z), family = gaussian(),
                            data = d3t))
  if (!inherits(f4t, "error")) {
    cat("tibble min", f4t$frame$linpreds[[1]]$mo[[1]]$min, "\n")
  }
  f4m <- try_("mo in sigma", frm(bf(y ~ z, sigma ~ mo(k)), family = gaussian(),
                                 data = d3))
  if (!inherits(f4m, "error")) {
    cat("sigma mo min", f4m$frame$linpreds[[2]]$mo[[1]]$min, "\n")
  }
  d3f <- d3
  attr(d3f$k, "min") <- 0.5
  try_("attr min 0.5", frm(bf(y ~ mo(k)), family = gaussian(), data = d3f))
  attr(d3f$k, "min") <- 2L
  try_("attr min 2", frm(bf(y ~ mo(k)), family = gaussian(), data = d3f))
}

cat("\n#### 2980 frm_multiple + int_conditions\n")
set.seed(12)
imps <- lapply(1:3, function(i) {
  dd <- rd[1:80, c("y", "x", "z")]
  dd$y <- dd$y + rnorm(80, sd = 0.01)
  dd
})
fm <- try_("frm_multiple", frm_multiple(bf(y ~ x * z), family = gaussian(),
                                         data = imps))
if (!inherits(fm, "error")) {
  ce <- try_("ce int_conditions", conditional_effects(
    fm, effects = "x:z", int_conditions = list(z = c(-1, 1))))
  if (!inherits(ce, "error")) {
    cat("z values:", unique(ce[[1]]$z), " rows", nrow(ce[[1]]), "\n")
  }
  ce2 <- try_("ce conditions", conditional_effects(
    fm, effects = "x", conditions = data.frame(z = c(-1, 1))))
  if (!inherits(ce2, "error")) {
    cat("conditions z:", unique(ce2[[1]]$z), " cond__",
        unique(as.character(ce2[[1]]$cond__)), "\n")
  }
  ce3 <- try_("ce int_conditions + boot arg", conditional_effects(
    fm, effects = "x:z", int_conditions = list(z = c(-1, 1)),
    resolution = 5))
  if (!inherits(ce3, "error")) cat("resolution 5 rows", nrow(ce3[[1]]), "\n")
}

cat("\n#### 4337 default_prior on rr()\n")
set.seed(13)
dr <- data.frame(g = factor(rep(1:10, each = 12)),
                 x1 = rnorm(120), x2 = rnorm(120), x3 = rnorm(120))
dr$y <- rnorm(120, 0.5 * dr$x1 + rnorm(10)[dr$g])
dp <- try_("default_prior rr", default_prior(
  bf(y ~ x1 + rr(0 + x1 + x2 + x3 | g, d = 1)), family = gaussian(),
  data = dr))
if (!inherits(dp, "error")) print(dp)
dp2 <- try_("default_prior rr + (1|g)", default_prior(
  bf(y ~ x1 + rr(0 + x1 + x2 + x3 | g, d = 1) + (1 | g)),
  family = gaussian(), data = dr))
if (!inherits(dp2, "error")) print(dp2[dp2$class %in% c("sd", "cor"), ])

cat("\n#### 5323 get_predicted on simplex families\n")
set.seed(5)
dm <- data.frame(x = rnorm(60))
pm <- cbind(1, exp(0.5 * dm$x), exp(-0.3 * dm$x))
dm$Y <- pm / rowSums(pm)
dm$Y <- dm$Y * 0.9 + 0.1 / 3
fd <- try_("dirichlet", frm(bf(Y ~ x), family = dirichlet(), data = dm))
if (!inherits(fd, "error")) {
  try_("get_predicted dirichlet", insight::get_predicted(fd))
  try_("get_predicted dirichlet link", insight::get_predicted(fd,
                                                              predict = "link"))
}
dm$yc <- factor(sample(c("a", "b", "c"), 60, TRUE))
fc <- try_("categorical", frm(bf(yc ~ x), family = categorical(), data = dm))
if (!inherits(fc, "error")) {
  r <- try_("get_predicted categorical", insight::get_predicted(fc))
  if (!inherits(r, "error")) cat("categorical: class", class(r)[1],
                                 " length", length(r), "\n")
}
