# Reviewer: offset() in emmeans() and conditional_effects() (claim 2).
# Data seed 21, the same construction as formrobust-rev-brms-emm.R.
# REVLIB="" for the base arm.
LIB <- Sys.getenv("REVLIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(emmeans)})
cat("frmtmb from", find.package("frmtmb"), "\n")
tr <- function(label, expr) {
  w <- character(0)
  r <- tryCatch(withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) cat(sprintf("[%s] ERROR: %s\n", label, r))
  if (length(w)) cat(sprintf("[%s] WARN: %s\n", label, unique(w)))
  invisible(if (inherits(r, "err")) NULL else r)
}
set.seed(21)
n <- 200
d <- data.frame(x = rnorm(n), time = runif(n, 1, 3),
                f = factor(sample(c("a", "b"), n, TRUE)),
                z = runif(n, 0, 1))
d$y <- rpois(n, exp(0.3 + 0.4 * d$x + 0.2 * (d$f == "b")) * d$time)
d$y2 <- rnorm(n, 1 + 0.5 * d$x + d$z, 0.5)
d$y3 <- rnorm(n, 0.5 * d$x, exp(0.2 + log(d$time)))
d$g <- factor(rep(1:20, each = 10))
mx <- mean(d$x); lmt <- log(mean(d$time))

fp <- frm(bf(y ~ x + f + offset(log(time))), data = d, family = poisson())
b <- fixef(fp)[, "Estimate"]
em <- tr("emmeans pois", summary(emmeans(fp, ~ f)))
if (!is.null(em)) {
  cat("emmean - (b0 + bx mean(x) [+ bfb]):",
      format(em$emmean - (b[["Intercept"]] + b[["x"]] * mx +
                            c(0, b[["fb"]])), digits = 4),
      " log(mean(time)) =", lmt, "\n")
}
em <- tr("emmeans pois epred", summary(emmeans(fp, ~ f, epred = TRUE)))
if (!is.null(em)) {
  cat("epred emmean / exp(lp + log(mean(time))):",
      format(em$emmean / exp(b[["Intercept"]] + b[["x"]] * mx +
                               c(0, b[["fb"]]) + lmt), digits = 15), "\n")
}
em <- tr("emmeans pois epred at time = 2",
         summary(emmeans(fp, ~ f, epred = TRUE, at = list(time = 2))))
if (!is.null(em)) {
  cat("epred at time=2 / exp(lp + log 2):",
      format(em$emmean / exp(b[["Intercept"]] + b[["x"]] * mx +
                               c(0, b[["fb"]]) + log(2)), digits = 15), "\n")
}
em <- tr("emmeans pois type=response",
         summary(emmeans(fp, ~ f, type = "response")))
if (!is.null(em)) print(em)
ce <- tr("CE pois", conditional_effects(fp, effects = "x"))
if (!is.null(ce)) {
  c1 <- ce[[1]]
  cat("CE names:", names(c1), "\n")
  cat("CE time column:", unique(c1$time), " mean(time)", mean(d$time), "\n")
  cat("CE / exp(b0 + bx x + log(mean time)) range:",
      format(range(c1$estimate__ / exp(b[["Intercept"]] + b[["x"]] * c1$x +
                                         lmt)), digits = 15), "\n")
}
ce <- tr("CE default effects", conditional_effects(fp))
if (!is.null(ce)) cat("CE default displays:", names(ce), "\n")
ce <- tr("CE with conditions time = 2",
         conditional_effects(fp, effects = "x",
                             conditions = data.frame(time = 2)))
if (!is.null(ce)) {
  c1 <- ce[[1]]
  cat("CE(time=2) / exp(b0+bx x+log 2) range:",
      format(range(c1$estimate__ / exp(b[["Intercept"]] + b[["x"]] * c1$x +
                                         log(2))), digits = 15), "\n")
}
ce <- tr("CE effects = time", conditional_effects(fp, effects = "time"))
if (!is.null(ce)) cat("CE time range", range(ce[[1]]$time), "\n")

cat("\n-- random effects (grid route?) --\n")
fr <- frm(bf(y ~ x + (1 | g) + offset(log(time))), data = d,
          family = poisson())
br <- fixef(fr)[, "Estimate"]
em <- tr("emmeans RE", summary(emmeans(fr, ~ 1)))
if (!is.null(em)) cat("RE emmean - (b0 + bx mx):",
                      em$emmean - (br[["Intercept"]] + br[["x"]] * mx), "\n")
em <- tr("emmeans RE epred", summary(emmeans(fr, ~ 1, epred = TRUE)))
if (!is.null(em)) cat("RE epred / exp(b0 + bx mx + lmt):",
                      format(em$emmean / exp(br[["Intercept"]] +
                                               br[["x"]] * mx + lmt),
                             digits = 15), "\n")

cat("\n-- nonlinear with offset in a nonlinear parameter --\n")
fnl <- tr("nl fit", frm(bf(y2 ~ a + b * x, a ~ 1 + offset(z), b ~ 1,
                           nl = TRUE), data = d))
if (!is.null(fnl)) {
  bn <- fixef(fnl)[, "Estimate"]
  print(bn)
  em <- tr("emmeans nl", summary(emmeans(fnl, ~ 1)))
  if (!is.null(em)) cat("nl emmean - (a + b mx):",
                        em$emmean - (bn[[1]] + bn[[2]] * mx),
                        " mean(z) =", mean(d$z), "\n")
  em <- tr("emmeans nl epred", summary(emmeans(fnl, ~ 1, epred = TRUE)))
  if (!is.null(em)) cat("nl epred - (a + b mx + mean z):",
                        em$emmean - (bn[[1]] + bn[[2]] * mx + mean(d$z)),
                        "\n")
  em <- tr("emmeans nlpar a", summary(emmeans(fnl, ~ 1, nlpar = "a")))
  if (!is.null(em)) cat("nlpar a emmean - a:", em$emmean - bn[[1]], "\n")
  ce <- tr("CE nl", conditional_effects(fnl, effects = "x"))
  if (!is.null(ce)) {
    c1 <- ce[[1]]
    cat("CE nl - (a + b x + mean z) range:",
        format(range(c1$estimate__ - (bn[[1]] + bn[[2]] * c1$x +
                                        mean(d$z))), digits = 4), "\n")
  }
}

cat("\n-- sigma offset --\n")
fs <- frm(bf(y3 ~ x, sigma ~ 1 + offset(log(time))), data = d)
bs <- fixef(fs)[, "Estimate"]
print(bs)
em <- tr("emmeans sigma", summary(emmeans(fs, ~ 1, dpar = "sigma")))
if (!is.null(em)) cat("sigma emmean - b_sigma:",
                      em$emmean - bs[["sigma_Intercept"]], "\n")
em <- tr("emmeans mu of sigma-offset model", summary(emmeans(fs, ~ 1)))
if (!is.null(em)) cat("mu emmean - (b0 + bx mx):",
                      em$emmean - (bs[["Intercept"]] + bs[["x"]] * mx), "\n")
ce <- tr("CE sigma model", conditional_effects(fs, effects = "x",
                                               dpar = "sigma"))
if (!is.null(ce)) cat("CE sigma / exp(b_s + log mean time) range:",
                      format(range(ce[[1]]$estimate__ /
                                     exp(bs[["sigma_Intercept"]] + lmt)),
                             digits = 15), "\n")

cat("\n-- offset with a factor inside, and offset(time) raw --\n")
d$ef <- factor(ifelse(d$time > 2, "hi", "lo"))
ff <- tr("offset(log(as.numeric(ef)))",
         frm(bf(y ~ x + offset(log(as.numeric(ef)))), data = d,
             family = poisson()))
if (!is.null(ff)) {
  tr("CE factor offset", conditional_effects(ff, effects = "x"))
  tr("emmeans factor offset", print(summary(emmeans(ff, ~ 1))))
  tr("emmeans factor offset epred", print(summary(emmeans(ff, ~ 1,
                                                          epred = TRUE))))
}
fo <- frm(bf(y ~ x + offset(time)), data = d, family = poisson())
tr("CE offset(time)", conditional_effects(fo, effects = "x"))
tr("emmeans offset(time)", print(summary(emmeans(fo, ~ 1))))
cat("names(fp$data):", names(fp$data), "\n")
tr("predict fp in-sample identical to newdata = d", {
  a <- fitted(fp); bb <- fitted(fp, newdata = d)
  cat("fitted in-sample vs newdata identical:", identical(unname(a[, 1]),
                                                          unname(bb[, 1])),
      "\n") })
# the offset variable enters na.action: NA in time drops the row
dn <- d; dn$time[3] <- NA
fn <- tr("offset var NA", frm(bf(y ~ x + offset(log(time))), data = dn,
                              family = poisson()))
if (!is.null(fn)) cat("nobs with time NA:", nobs(fn), "\n")
