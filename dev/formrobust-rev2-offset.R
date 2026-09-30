# Re-check (punch round 1): emmeans() offsets on shapes the lane did not
# try, against brms's emmeans() at FIXED parameters equal to frmtmb's
# estimates (the test helper brms_fixed_fit(), Fixed_param draws), and
# against hand formulas. Data seed 21. Stan cache: the reviewer's copy.
.libPaths(c("C:/Users/adf44/source/r/wt-formrobust-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/formrobust-rev-stan-cache",
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages({library(testthat); library(frmtmb); library(emmeans)})
tt <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/tests/testthat/"
for (h in c("helper-warnings.R", "helper-brms.R", "helper-brms-methods.R")) {
  sys.source(file.path(tt, h), envir = globalenv())
}
tr <- function(label, expr) {
  w <- character(0)
  r <- tryCatch(withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  }, message = function(c) invokeRestart("muffleMessage")),
  error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) cat(sprintf("[%s] ERROR: %s\n", label, r))
  if (length(w)) cat(sprintf("[%s] WARN: %s\n", label, unique(w)))
  invisible(if (inherits(r, "err")) NULL else r)
}
ulp <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300)) /
  .Machine$double.eps
cmp <- function(label, ef, eb) {
  if (is.null(ef) || is.null(eb)) return(invisible())
  a <- as.data.frame(ef); b <- as.data.frame(eb)
  cat(sprintf("[%s] frmtmb: %s\n", label,
              paste(format(a$emmean, digits = 10), collapse = " ")))
  cat(sprintf("[%s] brms  : %s\n", label,
              paste(format(b$emmean, digits = 10), collapse = " ")))
  if (nrow(a) == nrow(b)) {
    cat(sprintf("[%s] max rel diff in ulp: %.1f\n", label,
                ulp(a$emmean, b$emmean)))
  } else cat(sprintf("[%s] ROW COUNTS DIFFER %d %d\n", label, nrow(a),
                     nrow(b)))
}
set.seed(21)
n <- 200
d <- data.frame(x = rnorm(n), time = runif(n, 1, 3),
                f = factor(sample(c("a", "b"), n, TRUE)),
                h = factor(sample(c("u", "v"), n, TRUE)),
                z = runif(n, 0, 1))
d$yc <- rpois(n, exp(0.3 + 0.4 * d$x + 0.2 * (d$f == "b")) * d$time)
d$yc2 <- rpois(n, exp(0.1 + 0.3 * d$x + d$z) * d$time)
d$y3 <- rnorm(n, 0.5 * d$x + d$z, exp(0.2 + log(d$time)))
d$y2 <- rnorm(n, 1 + 0.5 * d$x)
mx <- mean(d$x); lmt <- log(mean(d$time))

cat("\n== A. offset + factors, by = and at = ==\n")
bfA <- bf(yc ~ x + f + h + offset(log(time)))
fA <- frm(bfA, data = d, family = poisson())
bA <- tr("brms fixed A", brms_fixed_fit(brms::bf(yc ~ x + f + h +
                                                   offset(log(time))),
                                        brms::brmsfamily("poisson"), d, fA,
                                        ndraws = 4))
for (sp in list(list("~ f by h", ~ f, list(by = "h")),
                list("~ f | h", ~ f | h, list()),
                list("~ f at time c(1,2)", ~ f,
                     list(at = list(time = c(1, 2)))),
                list("~ f | time at c(1,2)", ~ f | time,
                     list(at = list(time = c(1, 2)))),
                list("~ f epred at time c(1,2)", ~ f,
                     list(at = list(time = c(1, 2)), epred = TRUE)),
                list("~ x at x c(-1,1)", ~ x, list(at = list(x = c(-1, 1)))))) {
  ef <- tr(paste("frmtmb", sp[[1]]),
           do.call(emmeans, c(list(fA, sp[[2]]), sp[[3]])))
  eb <- if (!is.null(bA)) tr(paste("brms", sp[[1]]),
                             do.call(emmeans, c(list(bA, sp[[2]]), sp[[3]])))
  if (!is.null(ef) && !is.null(eb)) cmp(sp[[1]], summary(ef), summary(eb))
}
b <- fixef(fA)[, "Estimate"]
em <- summary(emmeans(fA, ~ f, at = list(time = c(1, 2))))
cat("hand: ~ f at time c(1,2) minus (lin + mean(log(c(1,2)))):",
    em$emmean - (b[["Intercept"]] + b[["x"]] * mx + c(0, b[["fb"]]) +
                   b[["hv"]] * mean(d$h == "v") * 0 + b[["hv"]] / 2 +
                   mean(log(c(1, 2)))), "\n")

cat("\n== B. two offsets in one formula ==\n")
fB <- frm(bf(yc2 ~ x + f + offset(log(time)) + offset(z)), data = d,
          family = poisson())
bB <- tr("brms fixed B", brms_fixed_fit(brms::bf(yc2 ~ x + f +
                                                   offset(log(time)) +
                                                   offset(z)),
                                        brms::brmsfamily("poisson"), d, fB,
                                        ndraws = 4))
for (ep in c(FALSE, TRUE)) {
  ef <- tr("frmtmb B", summary(emmeans(fB, ~ f, epred = ep)))
  eb <- if (!is.null(bB)) tr("brms B", summary(emmeans(bB, ~ f, epred = ep)))
  cmp(paste("two offsets, epred =", ep), ef, eb)
}
bb <- fixef(fB)[, "Estimate"]
ef <- summary(emmeans(fB, ~ f))
cat("hand: minus (lin + log(mean time) + mean z):",
    ef$emmean - (bb[["Intercept"]] + bb[["x"]] * mx + c(0, bb[["fb"]]) +
                   lmt + mean(d$z)), "\n")

cat("\n== C. offset in mu and in sigma ==\n")
bfC <- bf(y3 ~ x + offset(z), sigma ~ 1 + offset(log(time)))
fC <- frm(bfC, data = d)
bC <- tr("brms fixed C", brms_fixed_fit(brms::bf(y3 ~ x + offset(z),
                                                 sigma ~ 1 +
                                                   offset(log(time))),
                                        brms::brmsfamily("gaussian"), d, fC,
                                        ndraws = 4))
for (dp in list(NULL, "sigma")) {
  for (ep in c(FALSE, TRUE)) {
    if (ep && !is.null(dp)) next
    ef <- tr("frmtmb C", summary(emmeans(fC, ~ 1, dpar = dp, epred = ep)))
    eb <- if (!is.null(bC)) tr("brms C", summary(emmeans(bC, ~ 1, dpar = dp,
                                                         epred = ep)))
    cmp(paste("mu+sigma offsets, dpar =", dp %||% "mu", "epred =", ep),
        ef, eb)
  }
}
ef <- tr("frmtmb C at time", summary(emmeans(fC, ~ 1, dpar = "sigma",
                                             at = list(time = c(1, 2)))))
eb <- if (!is.null(bC)) tr("brms C at time",
                           summary(emmeans(bC, ~ 1, dpar = "sigma",
                                           at = list(time = c(1, 2)))))
cmp("sigma at time c(1,2)", ef, eb)

cat("\n== D. multivariate, offset on one response only ==\n")
bfD <- bf(yc ~ x + f + offset(log(time)), family = poisson()) +
  bf(y2 ~ x + f, family = gaussian())
fD <- frm(bfD, data = d)
bD <- tr("brms fixed D", brms_fixed_fit(
  brms::bf(yc ~ x + f + offset(log(time)), family = brms::brmsfamily("poisson")) +
    brms::bf(y2 ~ x + f, family = brms::brmsfamily("gaussian")) +
    brms::set_rescor(FALSE), NULL, d, fD, ndraws = 4))
bd <- fixef(fD)[, "Estimate"]
lin_yc <- bd[["yc_Intercept"]] + bd[["yc_x"]] * mx + c(0, bd[["yc_fb"]])
lin_y2 <- bd[["y2_Intercept"]] + bd[["y2_x"]] * mx + c(0, bd[["y2_fb"]])
specs <- list(list("both resp", list()),
              list("resp yc", list(resp = "yc")),
              list("resp y2", list(resp = "y2")),
              list("both, at time c(1,2)", list(at = list(time = c(1, 2)))),
              list("resp yc, at time 1", list(resp = "yc",
                                               at = list(time = 1))),
              list("both epred", list(epred = TRUE)))
for (sp in specs) {
  ef <- tr(paste("frmtmb D", sp[[1]]),
           summary(do.call(emmeans, c(list(fD, ~ f), sp[[2]]))))
  eb <- if (!is.null(bD)) tr(paste("brms D", sp[[1]]),
                             summary(do.call(emmeans,
                                             c(list(bD, ~ f), sp[[2]]))))
  cmp(paste("mv", sp[[1]]), ef, eb)
}
ef <- summary(emmeans(fD, ~ f))
cat("mv hand: frmtmb rows:", nrow(ef), "\n")
print(as.data.frame(ef)[, 1:3])
cat("hand yc: lin + lmt =", format(lin_yc + lmt, digits = 10),
    " y2: lin =", format(lin_y2, digits = 10), "\n")
ef <- summary(emmeans(fD, ~ f, at = list(time = c(1, 2))))
print(as.data.frame(ef)[, 1:3])
cat("hand yc at time c(1,2) averaged:", format(lin_yc + mean(log(c(1, 2))),
                                               digits = 10), "\n")
