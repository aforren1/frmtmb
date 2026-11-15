# Lane parity2: frmtmb on each item of the lane, on the data of
# dev/parity2-data.R. Arm "base" reads rellib-r21 (0.82.0) alone; arm
# "lane" puts the lane library in front.
#   Rscript dev/parity2-frm-probe.R base > dev/parity2-out/frm-base.txt 2>&1
#   Rscript dev/parity2-frm-probe.R lane > dev/parity2-out/frm-lane.txt 2>&1
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-parity2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, " lib:", find.package("frmtmb"),
    as.character(packageVersion("frmtmb")), "\n")
here <- "C:/Users/adf44/source/r/frmtmb-wt-parity2/dev"
source(file.path(here, "parity2-data.R"))
try_ <- function(label, expr) {
  w <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  if (inherits(r, "error")) cat(label, ": ERROR ", conditionMessage(r), "\n",
                                sep = "")
  for (x in w) cat("   warning: ", substr(x, 1, 200), "\n", sep = "")
  invisible(r)
}
ce_names <- function(label, fit, ...) {
  r <- try_(label, conditional_effects(fit, ...))
  if (!inherits(r, "error")) {
    cat(label, ": ", paste(names(r), collapse = " "), "\n", sep = "")
  }
  invisible(r)
}

cat("\n#### 3922 default effects\n")
shapes <- list(
  expxz = list(f = bf(y ~ exp(x):z), fam = gaussian()),
  t2 = list(f = bf(y ~ t2(x, z)), fam = gaussian()),
  moz = list(f = bf(y ~ mo(k):z), fam = gaussian()),
  cs = list(f = bf(yo ~ cs(x) + z), fam = sratio()),
  sig = list(f = bf(y ~ x, sigma ~ exp(z)), fam = gaussian()),
  mv = list(f = bf(y ~ x) + bf(y2 ~ z) + set_rescor(FALSE),
            fam = gaussian()),
  mi = list(f = bf(y ~ mi(x2) + z) + bf(x2 | mi() ~ w) + set_rescor(FALSE),
            fam = gaussian()),
  three = list(f = bf(y ~ x * z * w), fam = gaussian()),
  nl = list(f = bf(y ~ a * exp(b * x), a ~ 1 + w, b ~ z, nl = TRUE),
            fam = gaussian()),
  re = list(f = bf(y ~ x + (1 + z | g), sigma ~ w), fam = gaussian()),
  polyk = list(f = bf(y ~ poly(x, 2) + z + x:w), fam = gaussian())
)
for (nm in names(shapes)) {
  s <- shapes[[nm]]
  fit <- try_(nm, frm(s$f, family = s$fam, data = d))
  if (!inherits(fit, "error")) ce_names(paste(nm, "ce"), fit)
}

cat("\n#### 3919 s(exp(x))\n")
fs <- try_("sexp", frm(bf(y ~ s(exp(x))), family = gaussian(), data = d))
if (!inherits(fs, "error")) {
  cat("variables:", paste(variables(fs), collapse = " "), "\n")
  ce_names("sexp ce", fs)
  cat("fitted at x = 1:", fitted(fs, newdata = data.frame(x = 1))[1, 1],
      "\n")
}

cat("\n#### 4160 attr(x, 'min') of a mo() variable\n")
set.seed(9)
d3 <- data.frame(k = rep(1:4, 10))
d3$y <- 0.5 * d3$k + rnorm(40)
attr(d3$k, "min") <- 0
f3 <- try_("mo attr min 0", frm(bf(y ~ mo(k)), family = gaussian(),
                                 data = d3))
if (!inherits(f3, "error")) {
  m <- f3$frame$linpreds[[1]]$mo[[1]]
  cat("frame mo min", m$min, " D", m$D, "\n")
  cat("fitted at k = 1..4:", format(fitted(f3, newdata = data.frame(k = 1:4))[,
      1], digits = 6), "\n")
  cat("fitted in-sample rows 1..4:", format(fitted(f3)[1:4, 1], digits = 6),
      "\n")
}

cat("\n#### 1560 lf() on a dpar bf() already sets\n")
f15 <- try_("lf", bf(y ~ x, sigma ~ z) + lf(sigma ~ w))
if (!inherits(f15, "error")) print(f15)

cat("\n#### 3373 REML AIC against lme4\n")
data("sleepstudy", package = "lme4")
fr <- frm(bf(Reaction ~ Days + (Days | Subject)), family = gaussian(),
          data = sleepstudy, REML = TRUE)
lr <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = TRUE)
cat("frmtmb AIC", AIC(fr), " df", attr(logLik(fr), "df"), " logLik",
    as.numeric(logLik(fr)), "\n")
cat("lme4   AIC", AIC(lr), " df", attr(logLik(lr), "df"), " logLik",
    as.numeric(logLik(lr)), "\n")
cat("frmtmb BIC", BIC(fr), " lme4 BIC", BIC(lr), "\n")

cat("\n#### 5337 drop_unused_levels = FALSE\n")
fdu <- try_("dul", frm(bf(y ~ f + (1 | g)), family = gaussian(), data = du,
                       drop_unused_levels = FALSE))
if (!inherits(fdu, "error")) {
  cat("fixef:", names(fixef(fdu)[, 1]), "\n")
  cat("ranef rows:", rownames(ranef(fdu)$g), "\n")
  cat("variables:", paste(grep("^(b_|r_g|sd_)", variables(fdu), value = TRUE),
                          collapse = " "), "\n")
}

cat("\n#### 5315 b_Intercept order under 0 + Intercept\n")
fi <- try_("icpt", frm(bf(y ~ 0 + Intercept + x, sigma ~ z),
                       family = gaussian(), data = d))
if (!inherits(fi, "error")) cat("variables:", paste(variables(fi),
                                                   collapse = " "), "\n")
fi2 <- try_("icpt mu only", frm(bf(y ~ 0 + Intercept + x),
                                family = gaussian(), data = d))
if (!inherits(fi2, "error")) cat("variables:", paste(variables(fi2),
                                                    collapse = " "), "\n")
fi3 <- try_("icpt sigma 0+I", frm(bf(y ~ x, sigma ~ 0 + Intercept + z),
                                  family = gaussian(), data = d))
if (!inherits(fi3, "error")) cat("variables:", paste(variables(fi3),
                                                    collapse = " "), "\n")

cat("\n#### 5323 insight::get_predicted on multinomial\n")
fmn <- try_("multinomial", frm(bf(Y | trials(size) ~ x),
                               family = multinomial(), data = dm))
if (!inherits(fmn, "error")) {
  gp <- try_("get_predicted", insight::get_predicted(fmn))
  if (!inherits(gp, "error")) {
    cat("class", class(gp), " length", length(gp), "\n")
  }
}

cat("\n#### 2980 frm_multiple + int_conditions\n")
fm <- try_("frm_multiple", frm_multiple(bf(y ~ x * z), family = gaussian(),
                                        data = list(d, d)))
if (!inherits(fm, "error")) {
  r <- ce_names("multiple ce int_conditions", fm, effects = "x:z",
                int_conditions = list(z = c(-1, 1)))
  if (!inherits(r, "error")) print(unique(r[[1]]$z))
}

cat("\n#### 4337 default_prior on rr()\n")
set.seed(2)
drr <- data.frame(g = factor(rep(1:10, each = 6)),
                  tf = factor(rep(c("a", "b", "c"), 20)), x = rnorm(60))
drr$y <- rnorm(60)
dp4 <- try_("default_prior rr", default_prior(
  bf(y ~ x + rr(0 + tf | g, d = 1)), family = gaussian(), data = drr))
if (!inherits(dp4, "error")) print(as.data.frame(dp4)[, 1:4])

cat("\n#### 1497 new levels: y ~ x + (1 | g), band = boot\n")
fg <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dg)
ceb <- try_("boot", conditional_effects(fg, "x", re_formula = NULL,
                                        band = "boot", boot = 200, seed = 2,
                                        int_conditions = list(x = c(-1, 0,
                                                                    1))))
if (!inherits(ceb, "error")) {
  cat("boot width at x = -1, 0, 1:",
      format(ceb$x$upper__ - ceb$x$lower__, digits = 4), "\n")
}
cew <- try_("wald", conditional_effects(fg, "x", re_formula = NULL,
                                        int_conditions = list(x = c(-1, 0,
                                                                    1))))
if (!inherits(cew, "error")) {
  cat("wald width at x = -1, 0, 1:",
      format(cew$x$upper__ - cew$x$lower__, digits = 4), "\n")
}
cat("done\n")
