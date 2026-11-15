# Review of lane parity2, item 3373: the df of a REML fit. lme4 equality
# on more models, the readers of df, and the invariant df(REML) ==
# df(ML) that lme4 keeps, on shapes where the beta template could hold
# columns the fit does not estimate (an aliased column, a zero column
# of an unused level, a mo() term on one value, a sigma formula).
#   Rscript dev/reviews/parity2-review-reml.R base|lane
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-parity2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(lme4)})
cat("arm", arm, " lib:", find.package("frmtmb"),
    as.character(packageVersion("frmtmb")), "\n")
q <- function(expr) suppressMessages(suppressWarnings(expr))
cmp <- function(label, ff, fl) {
  lf <- logLik(ff); ll <- logLik(fl)
  cat(sprintf(
    "%-12s frmtmb logLik %.4f df %s AIC %.3f BIC %.3f nobs %s dfres %s | lme4 logLik %.4f df %s AIC %.3f BIC %.3f nobs %s dfres %s\n",
    label, as.numeric(lf), attr(lf, "df"), AIC(ff), BIC(ff), nobs(ff),
    df.residual(ff), as.numeric(ll), attr(ll, "df"), AIC(fl), BIC(fl),
    nobs(fl), df.residual(fl)))
  ea <- stats::extractAIC(ff)
  cat(sprintf("%-12s extractAIC frmtmb %s %.3f\n", label, ea[1], ea[2]))
}
data("sleepstudy", package = "lme4")
data("Penicillin", package = "lme4")
data("cake", package = "lme4")

cat("\n#### lme4 equality\n")
cmp("rslope", q(frm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
                   REML = TRUE)),
    lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML = TRUE))
cmp("uncorr", q(frm(Reaction ~ Days + (Days || Subject), data = sleepstudy,
                   REML = TRUE)),
    lmer(Reaction ~ Days + (Days || Subject), sleepstudy, REML = TRUE))
cmp("twofac", q(frm(diameter ~ 1 + (1 | plate) + (1 | sample),
                   data = Penicillin, REML = TRUE)),
    lmer(diameter ~ 1 + (1 | plate) + (1 | sample), Penicillin, REML = TRUE))
cmp("cake18", q(frm(angle ~ recipe * temperature + (1 | recipe:replicate),
                   data = cake, REML = TRUE)),
    lmer(angle ~ recipe * temperature + (1 | recipe:replicate), cake,
         REML = TRUE))
ss <- sleepstudy
ss$Days2 <- 2 * ss$Days
cmp("aliased", q(frm(Reaction ~ Days + Days2 + (1 | Subject), data = ss,
                    REML = TRUE)),
    q(lmer(Reaction ~ Days + Days2 + (1 | Subject), ss, REML = TRUE)))

cat("\n#### anova of two REML fits, same fixed effects, against lme4\n")
r0 <- q(frm(Reaction ~ Days + (1 | Subject), data = sleepstudy, REML = TRUE))
r1 <- q(frm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
            REML = TRUE))
print(anova(r0, r1))
print(anova(lmer(Reaction ~ Days + (1 | Subject), sleepstudy),
            lmer(Reaction ~ Days + (Days | Subject), sleepstudy),
            refit = FALSE))

cat("\n#### weights: df only (lme4's weights are precision weights)\n")
set.seed(3)
ss$wt <- sample(1:3, nrow(ss), TRUE)
fw <- q(frm(Reaction | weights(wt) ~ Days + (1 | Subject), data = ss,
            REML = TRUE))
fwm <- q(frm(Reaction | weights(wt) ~ Days + (1 | Subject), data = ss,
             REML = FALSE))
cat("weights REML df", attr(logLik(fw), "df"), " ML df",
    attr(logLik(fwm), "df"), " lme4 REML df",
    attr(logLik(lmer(Reaction ~ Days + (1 | Subject), ss, weights = wt)),
         "df"), "\n")

cat("\n#### df(REML) == df(ML) on shapes with a non-estimated beta\n")
inv <- function(label, f, data, ...) {
  a <- tryCatch(q(frm(f, data = data, REML = TRUE, ...)),
                error = function(e) e)
  b <- tryCatch(q(frm(f, data = data, REML = FALSE, ...)),
                error = function(e) e)
  if (inherits(a, "error") || inherits(b, "error")) {
    cat(label, "ERROR", conditionMessage(if (inherits(a, "error")) a else b),
        "\n")
    return(invisible())
  }
  da <- tryCatch(attr(logLik(a), "df"), error = function(e) conditionMessage(e))
  db <- tryCatch(attr(logLik(b), "df"), error = function(e) conditionMessage(e))
  cat(sprintf("%-10s REML df %s  ML df %s  beta len %d  opt REML %d ML %d\n",
              label, da, db, length(a$frame$par_template$beta),
              length(a$opt$par), length(b$opt$par)))
}
set.seed(4)
du <- data.frame(g = factor(sample(1:6, 60, TRUE)),
                 f = factor(sample(c("a", "b"), 60, TRUE),
                            levels = c("a", "b", "c")),
                 x = rnorm(60), z = rnorm(60), k = 2L)
du$k[1:30] <- rep(1:3, 10)
du$y <- rnorm(60, 1 + 0.4 * du$x + rnorm(6)[du$g])
inv("plain", y ~ x + (1 | g), du)
inv("aliased", y ~ x + I(2 * x) + (1 | g), du)
inv("unused", y ~ x + f + (1 | g), du, drop_unused_levels = FALSE)
inv("sigma", bf(y ~ x + (1 | g), sigma ~ z), du)
du1 <- du
du1$k <- 2L
inv("mo1val", y ~ x + mo(k) + (1 | g), du1)
inv("mo", y ~ x + mo(k) + (1 | g), du)
