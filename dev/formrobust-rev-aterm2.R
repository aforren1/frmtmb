# Reviewer: addition-term expressions, second pass. Seed 101 (same data
# as formrobust-rev-aterm.R). REVLIB="" for the base arm.
LIB <- Sys.getenv("REVLIB", "C:/Users/adf44/source/r/wt-formrobust-lib")
.libPaths(c(if (nzchar(LIB)) LIB, "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
bf <- frmtmb::bf
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
same <- function(label, a, b) {
  cat(sprintf("[%s] identical: %s\n", label, identical(a, b)))
}
sd_ <- function(...) suppressMessages(suppressWarnings(brms::standata(...)))
set.seed(101)
n <- 60
d <- data.frame(x = rnorm(n), t = runif(n, 0.5, 2), c = rpois(n, 4) + 1L,
                wt = runif(n, 0.5, 2), s = runif(n, 0.2, 0.6),
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                s1 = rep(c(TRUE, FALSE), length.out = n))
d$y <- rnorm(n, 0.5 * d$x)
d$y2 <- rnorm(n, -0.5 * d$x)
d$yc <- rpois(n, exp(0.2 + 0.3 * d$x) * d$t)
d$yb <- rbinom(n, d$c + 1, plogis(0.3 * d$x))

cat("\n== A. readers of the frame (weights(wt*k), env k) ==\n")
k <- 3
d$wk <- d$wt * 3
fa <- frm(bf(y | weights(wt * k) ~ x + (1 | f)), data = d)
fb <- frm(bf(y | weights(wk) ~ x + (1 | f)), data = d)
same("logLik", logLik(fa), logLik(fb))
tr("residuals", same("residuals", residuals(fa), residuals(fb)))
tr("influence", same("influence", influence(fa, groups = "f"),
                     influence(fb, groups = "f")))
tr("bootstrap", same("bootstrap", frm_bootstrap(fa, nsim = 3, seed = 1),
                     frm_bootstrap(fb, nsim = 3, seed = 1)))
tr("refit", same("refit", logLik(refit(fa, simulate(fa, seed = 1)[[1]])),
                 logLik(refit(fb, simulate(fb, seed = 1)[[1]]))))
tr("predict newdata", {
  set.seed(5); pa <- predict(fa, newdata = d[1:5, ])
  set.seed(5); pb <- predict(fb, newdata = d[1:5, ])
  same("predict newdata", pa, pb) })
tr("conditional_effects", same("ce", conditional_effects(fa)[[1]]$estimate__,
                               conditional_effects(fb)[[1]]$estimate__))
# k removed from the calling scope after the fit: the formula env holds it
rm(k)
tr("influence after rm(k)", {
  x1 <- influence(fa, groups = "f"); cat("  ok\n") })

cat("\n== B. trunc(min(y)) refits on subsets: brms re-evaluates ==\n")
d$lbm <- min(d$y) - 1
fa <- frm(bf(y | trunc(lb = min(y) - 1) ~ x + (1 | f)), data = d)
fb <- frm(bf(y | trunc(lb = lbm) ~ x + (1 | f)), data = d)
same("logLik", logLik(fa), logLik(fb))
tr("bootstrap trunc(min(y))", {
  ba <- frm_bootstrap(fa, nsim = 3, seed = 1)
  bb <- frm_bootstrap(fb, nsim = 3, seed = 1)
  same("bootstrap trunc(min(y)) vs fixed lb", ba, bb)
  print(ba$t); print(bb$t) })

cat("\n== C. cens with a single-value or literal interval bound ==\n")
dc <- d
dc$cc <- sample(c("none", "interval", "right"), n, TRUE)
dc$lo <- dc$y - 0.5
dc$y_c <- ifelse(dc$cc == "interval", dc$lo, dc$y)
dc$ubc <- max(dc$y_c) + 10
fa <- tr("cens(cc, max(y_c) + 10)",
         frm(bf(y_c | cens(cc, max(y_c) + 10) ~ x), data = dc))
fl <- tr("cens(cc, 10)", frm(bf(y_c | cens(cc, 10) ~ x), data = dc))
fb <- tr("cens(cc, ubc)", frm(bf(y_c | cens(cc, ubc) ~ x), data = dc))
s <- tr("brms cens(cc, max(y_c) + 10)",
        sd_(y_c | cens(cc, max(y_c) + 10) ~ x, data = dc))
if (!is.null(s)) cat("  brms rcens head", head(s$rcens), "\n")
s <- tr("brms cens(cc, 10)", sd_(y_c | cens(cc, 10) ~ x, data = dc))
if (!is.null(s)) cat("  brms rcens head", head(s$rcens), "\n")
s2 <- tr("brms cens(cc, ubc)", sd_(y_c | cens(cc, ubc) ~ x, data = dc))
if (!is.null(s) && !is.null(s2)) same("brms scalar vs column", s$rcens,
                                      s2$rcens)
if (!is.null(fl)) same("cens literal vs column", logLik(fl), logLik(fb))

cat("\n== D. NA-yielding expressions: base drops the row ==\n")
for (lab in c("weights", "se", "trunc", "cens", "rate", "trials")) {
  f <- switch(lab,
    weights = y | weights(ifelse(x > 1, NA, wt)) ~ x,
    se = y | se(ifelse(x > 1, NA, s)) ~ x,
    trunc = y | trunc(lb = ifelse(x > 1, NA, -5)) ~ x,
    cens = y | cens(ifelse(x > 1, NA, "none")) ~ x,
    rate = yc | rate(ifelse(x > 1, NA, t)) ~ x,
    trials = yb | trials(ifelse(x > 1, NA, c)) ~ x)
  fam <- switch(lab, rate = poisson(), trials = binomial(), gaussian())
  r <- tr(paste("NA", lab), frm(bf(f), data = d, family = fam))
  if (!is.null(r)) cat(sprintf("  [%s] nobs %d logLik %s\n", lab, nobs(r),
                               format(logLik(r))))
}
cat("rows with x > 1:", sum(d$x > 1), "\n")

cat("\n== E. index() and mi(idx) with an expression ==\n")
dm <- data.frame(id = 1:30, x = rnorm(30))
dm$z <- dm$x + rnorm(30, sd = 0.3)
dm$z[c(4, 9)] <- NA
dm$y <- 1 + 0.5 * dm$x + rnorm(30)
dm$id3 <- dm$id * 3
fa <- tr("index(id*3)", frm(bf(y ~ mi(z, idx = id3)) +
                              bf(z | mi() + index(id * 3) ~ x), data = dm))
fb <- tr("index(id3)", frm(bf(y ~ mi(z, idx = id3)) +
                             bf(z | mi() + index(id3) ~ x), data = dm))
if (!is.null(fa)) same("index expr logLik", logLik(fa), logLik(fb))

cat("\n== F. weights(scale = TRUE) with NA rows and on newdata ==\n")
dn <- d; dn$x[1:4] <- NA
fa <- frm(bf(y | weights(wt, scale = TRUE) ~ x), data = dn)
s <- sd_(y | weights(wt, scale = TRUE) ~ x, data = dn)
same("scale=TRUE with NA rows vs brms", fa$frame$aterm_values[[1]]$weights,
     as.numeric(s$weights))
cat("  sum of weights", sum(fa$frame$aterm_values[[1]]$weights), "n",
    nobs(fa), "\n")
tr("scale=TRUE formula print", print(formula(fa)))
tr("scale=TRUE summary", invisible(capture.output(summary(fa))))
tr("update_adterms", print(update_adterms(bf(y | weights(wt, scale = TRUE) ~
                                                 x), ~ se(s))$formula))
tr("weights(scale=TRUE) positional", {
  f2 <- frm(bf(y | weights(wt, TRUE) ~ x), data = dn)
  same("positional scale", logLik(f2), logLik(fa)) })
tr("weights(x = wt)", frm(bf(y | weights(x = wt) ~ x), data = dn))
tr("weights()", frm(bf(y | weights() ~ x), data = dn))
tr("weights(wt, sc = TRUE)", frm(bf(y | weights(wt, sc = TRUE) ~ x),
                                 data = dn))
tr("brms weights(wt, sc = TRUE)", sd_(y | weights(wt, sc = TRUE) ~ x,
                                      data = dn))
tr("weights(wt, FALSE, 3)", frm(bf(y | weights(wt, FALSE, 3) ~ x),
                                data = dn))
