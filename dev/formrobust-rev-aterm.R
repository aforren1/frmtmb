# Reviewer: attacks on addition-term expressions (lane claim 1).
# Seed 101. REVLIB="" for the base arm (rellib-r4).
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
  if (length(w)) cat(sprintf("[%s] WARN: %s\n", label, w))
  invisible(if (inherits(r, "err")) NULL else r)
}
same <- function(label, a, b) {
  cat(sprintf("[%s] identical: %s\n", label, identical(a, b)))
  if (!identical(a, b) && is.numeric(a) && is.numeric(b) &&
      length(a) == length(b)) {
    cat(sprintf("    max rel diff %.3g\n",
                max(abs(a - b) / pmax(abs(b), 1e-300))))
  }
}
sd_ <- function(...) {
  suppressMessages(suppressWarnings(brms::standata(...)))
}
fv <- function(fit) fit$frame$aterm_values[[1L]]

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

cat("\n== 1. names that are functions: t and c ==\n")
d$t2 <- d$t * 2; d$c1 <- d$c + 1; d$c2 <- d$c / 2
fa <- tr("weights(t*2)", frm(bf(y | weights(t * 2) ~ x), data = d))
fb <- tr("weights(t2)", frm(bf(y | weights(t2) ~ x), data = d))
if (!is.null(fa)) same("weights(t*2) logLik", logLik(fa), logLik(fb))
fa <- tr("trials(c+1)", frm(bf(yb | trials(c + 1) ~ x), data = d,
                            family = binomial()))
fb <- tr("trials(c1)", frm(bf(yb | trials(c1) ~ x), data = d,
                           family = binomial()))
if (!is.null(fa)) same("trials(c+1) logLik", logLik(fa), logLik(fb))
fa <- tr("se(c/2)", frm(bf(y | se(c / 2) ~ x), data = d))
fb <- tr("se(c2)", frm(bf(y | se(c2) ~ x), data = d))
if (!is.null(fa)) same("se(c/2) logLik", logLik(fa), logLik(fb))
fa <- tr("rate(t*2)", frm(bf(yc | rate(t * 2) ~ x), data = d,
                          family = poisson()))
if (!is.null(fa)) {
  s <- sd_(yc | rate(t * 2) ~ x, data = d, family = stats::poisson())
  same("rate(t*2) vs brms denom", fv(fa)[["rate"]], as.numeric(s$denom))
}
# a variable named like a function that is NOT in data: t() is base::t
dd <- d; dd$t <- NULL
tr("weights(t * 2), no column t", frm(bf(y | weights(t * 2) ~ x), data = dd))

cat("\n== 2. formula environment ==\n")
f_env <- local({
  k <- 3
  wv <- runif(n, 1, 2)
  list(k = k, wv = wv,
       f1 = y | weights(wt * k) ~ x,
       f2 = y | weights(wv * 2) ~ x)
})
d$wk <- d$wt * 3
fa <- tr("weights(wt*k) local k", frm(bf(f_env$f1), data = d))
fb <- tr("weights(wk)", frm(bf(y | weights(wk) ~ x), data = d))
if (!is.null(fa)) same("wt*k logLik", logLik(fa), logLik(fb))
# an environment vector with rows dropped by NA in x
dna <- d; dna$x[c(3, 7)] <- NA
d_wv <- dna; d_wv$wv2 <- f_env$wv * 2
fa <- tr("weights(wv*2) env vector, NA rows", frm(bf(f_env$f2), data = dna))
fb <- tr("weights(wv2) column, NA rows", frm(bf(y | weights(wv2) ~ x),
                                               data = d_wv))
if (!is.null(fa)) {
  same("env-vector weights logLik", logLik(fa), logLik(fb))
  same("env-vector weights values", fv(fa)[["weights"]],
       fv(fb)[["weights"]])
  s <- tr("brms env vector", sd_(f_env$f2, data = dna))
  if (!is.null(s)) same("env-vector vs brms", fv(fa)[["weights"]],
                        as.numeric(s$weights))
}
# influence / refit / bootstrap on the env-constant model
if (!is.null(fa <- tr("wt*k", frm(bf(f_env$f1), data = d)))) {
  tr("influence wt*k", {
    inf <- influence(fa, groups = "f"); cat("  influence ok\n") })
  tr("residuals wt*k", {
    r1 <- residuals(fa); r2 <- residuals(fb)
    same("residuals wt*k vs wk", r1, r2) })
  tr("bootstrap wt*k", {
    b1 <- frm_bootstrap(fa, nboot = 3, seed = 1)
    b2 <- frm_bootstrap(fb, nboot = 3, seed = 1)
    same("bootstrap wt*k vs wk", b1$t, b2$t) })
  tr("loo wt*k", { l <- loo(fa); cat("  loo elpd", l$estimates[1, 1], "\n") })
}

cat("\n== 3. factors and NA-yielding expressions ==\n")
d$fw <- as.numeric(d$f)
fa <- tr("weights(as.numeric(f))", frm(bf(y | weights(as.numeric(f)) ~ x),
                                        data = d))
fb <- tr("weights(fw)", frm(bf(y | weights(fw) ~ x), data = d))
if (!is.null(fa)) same("factor expr logLik", logLik(fa), logLik(fb))
fa <- tr("trials(as.integer(f) + c)", frm(bf(yb | trials(as.integer(f) + c) ~ x),
                                          data = d, family = binomial()))
fa <- tr("weights(f == 'a') logical", frm(bf(y | weights(f == "a") ~ x),
                                          data = d))
if (!is.null(fa)) {
  s <- tr("brms weights(f == 'a')", sd_(y | weights(f == "a") ~ x, data = d))
  if (!is.null(s)) same("weights(f=='a') vs brms", fv(fa)[["weights"]],
                        as.numeric(s$weights))
}
fa <- tr("weights(ifelse(x > 1, NA, wt))",
         frm(bf(y | weights(ifelse(x > 1, NA, wt)) ~ x), data = d))
if (!is.null(fa)) cat("  logLik", format(logLik(fa)), "\n")
tr("brms weights NA", sd_(y | weights(ifelse(x > 1, NA, wt)) ~ x, data = d))
fa <- tr("se(ifelse(x > 1, NA, s))",
         frm(bf(y | se(ifelse(x > 1, NA, s)) ~ x), data = d))
if (!is.null(fa)) cat("  logLik", format(logLik(fa)), "\n")
tr("brms se NA", sd_(y | se(ifelse(x > 1, NA, s)) ~ x, data = d))
fa <- tr("trunc(lb = ifelse(x > 1, NA, -5))",
         frm(bf(y | trunc(lb = ifelse(x > 1, NA, -5)) ~ x), data = d))
if (!is.null(fa)) cat("  logLik", format(logLik(fa)), "\n")
tr("brms trunc NA", sd_(y | trunc(lb = ifelse(x > 1, NA, -5)) ~ x, data = d))
fa <- tr("rate(ifelse(x > 1, NA, t))",
         frm(bf(yc | rate(ifelse(x > 1, NA, t)) ~ x), data = d,
             family = poisson()))
fa <- tr("trials(ifelse(x > 1, NA, c))",
         frm(bf(yb | trials(ifelse(x > 1, NA, c)) ~ x), data = d,
             family = binomial()))
if (!is.null(fa)) cat("  logLik", format(logLik(fa)), "\n")
# weights with an expression variable NA: row dropped (brms allvars)
dwn <- d; dwn$wt[5] <- NA
fa <- tr("weights(wt*2), wt NA in row 5",
         frm(bf(y | weights(wt * 2) ~ x), data = dwn))
if (!is.null(fa)) cat("  nobs", nobs(fa), "\n")
s <- tr("brms weights(wt*2) wt NA", sd_(y | weights(wt * 2) ~ x, data = dwn))
if (!is.null(s)) cat("  brms N", s$N, "\n")

cat("\n== 4. multivariate with subset() on one response ==\n")
fa <- tr("mv trunc(min(y)) + subset(s1) | y2 weights(wt/sum(wt))",
         frm(bf(y | trunc(lb = min(y) - 1) + subset(s1) ~ x) +
               bf(y2 | weights(wt / sum(wt)) ~ x), data = d))
s <- tr("brms mv", sd_(brms::bf(y | trunc(lb = min(y) - 1) + subset(s1) ~ x) +
                         brms::bf(y2 | weights(wt / sum(wt)) ~ x) +
                         brms::set_rescor(FALSE), data = d))
if (!is.null(fa) && !is.null(s)) {
  same("mv lb_y vs brms", fa$frame$aterm_values[["y"]][["trunc_lb"]],
       as.numeric(s$lb_y))
  same("mv weights_y2 vs brms", fa$frame$aterm_values[["y2"]][["weights"]],
       as.numeric(s$weights_y2))
}
fa <- tr("mv weights(wt/sum(wt)) + subset(s1) on y",
         frm(bf(y | weights(wt / sum(wt)) + subset(s1) ~ x) +
               bf(y2 | se(s / 2) ~ x), data = d))
s <- tr("brms mv2", sd_(brms::bf(y | weights(wt / sum(wt)) + subset(s1) ~ x) +
                          brms::bf(y2 | se(s / 2) ~ x) +
                          brms::set_rescor(FALSE), data = d))
if (!is.null(fa) && !is.null(s)) {
  same("mv2 weights_y vs brms", fa$frame$aterm_values[["y"]][["weights"]],
       as.numeric(s$weights_y))
  same("mv2 se_y2 vs brms", fa$frame$aterm_values[["y2"]][["se"]],
       as.numeric(s$se_y2))
}
fa <- tr("mv weights(wt, scale=TRUE) + subset(s1)",
         frm(bf(y | weights(wt, scale = TRUE) + subset(s1) ~ x) +
               bf(y2 ~ x), data = d))
s <- tr("brms mv3", sd_(brms::bf(y | weights(wt, scale = TRUE) +
                                   subset(s1) ~ x) +
                          brms::bf(y2 ~ x) + brms::set_rescor(FALSE),
                        data = d))
if (!is.null(fa) && !is.null(s)) {
  same("mv3 weights_y vs brms", fa$frame$aterm_values[["y"]][["weights"]],
       as.numeric(s$weights_y))
}
# the subset response's expression var with NA on a row it does not use
dsn <- d; dsn$wt[!dsn$s1][1:3] <- NA
fa <- tr("mv NA in wt outside the subset",
         frm(bf(y | weights(wt * 2) + subset(s1) ~ x) + bf(y2 ~ x),
             data = dsn))
s <- tr("brms mv NA outside subset",
        sd_(brms::bf(y | weights(wt * 2) + subset(s1) ~ x) +
              brms::bf(y2 ~ x) + brms::set_rescor(FALSE), data = dsn))
if (!is.null(fa) && !is.null(s)) {
  cat("  frmtmb N_y", length(fa$frame$aterm_values[["y"]][["weights"]]),
      " N_y2 ", length(fa$frame$y[["y2"]]), "; brms N_y", s$N_y, " N_y2",
      s$N_y2, "\n")
}

cat("\n== 5. newdata without an expression variable ==\n")
fr <- frm(bf(yc | rate(t * 2) ~ x), data = d, family = poisson())
nd <- d[1:5, c("x", "yc")]
tr("fitted rate(t*2), newdata lacks t", print(fitted(fr, newdata = nd)[, 1]))
tr("predict rate(t*2), newdata lacks t", predict(fr, newdata = nd))
tr("log_lik rate(t*2), newdata lacks t", log_lik(fr, newdata = nd))
fw <- frm(bf(y | weights(wt * 2) ~ x), data = d)
tr("fitted weights(wt*2), newdata lacks wt",
   print(fitted(fw, newdata = d[1:5, c("x", "y")])[, 1]))
tr("log_lik weights(wt*2), newdata lacks wt",
   print(log_lik(fw, newdata = d[1:5, c("x", "y")])[1, ]))
fs <- frm(bf(y | se(s / 2) ~ x), data = d)
tr("predict se(s/2), newdata lacks s",
   predict(fs, newdata = d[1:5, c("x", "y")]))
# the base arm's error for the bare-name spelling, for comparison
d$tt2 <- d$t * 2
fr2 <- frm(bf(yc | rate(tt2) ~ x), data = d, family = poisson())
tr("fitted rate(tt2), newdata lacks tt2", fitted(fr2, newdata = nd))

cat("\n== 6. cens() with an interval expression ==\n")
dc <- d
dc$lo <- dc$y - 0.5
dc$cc <- sample(c("none", "interval", "right"), n, TRUE)
dc$hi <- ifelse(dc$cc == "interval", dc$y + 0.5, NA)
dc$hi2 <- dc$hi * 1.5
dc$y_c <- ifelse(dc$cc == "interval", dc$lo, dc$y)
fa <- tr("cens(cc, hi * 1.5)", frm(bf(y_c | cens(cc, hi * 1.5) ~ x),
                                   data = dc))
fb <- tr("cens(cc, hi2)", frm(bf(y_c | cens(cc, hi2) ~ x), data = dc))
if (!is.null(fa)) same("cens interval expr logLik", logLik(fa), logLik(fb))
s <- tr("brms cens(cc, hi*1.5)", sd_(y_c | cens(cc, hi * 1.5) ~ x, data = dc))
if (!is.null(fa) && !is.null(s)) {
  cat("  brms rcens/Jint names:", paste(names(s)[grepl("cens|rcens|Jcens|y2|Y",
                                                       names(s))],
                                        collapse = " "), "\n")
}
# a single-value interval bound: brms recycles it
dc$y_c2 <- ifelse(dc$cc == "interval", dc$lo, dc$y)
fa <- tr("cens(cc, max(y) + 10)", frm(bf(y_c | cens(cc, max(y_c) + 10) ~ x),
                                      data = dc))
dc$ub_const <- max(dc$y_c) + 10
fb <- tr("cens(cc, ub_const)", frm(bf(y_c | cens(cc, ub_const) ~ x),
                                   data = dc))
if (!is.null(fa) && !is.null(fb)) {
  same("cens scalar interval bound logLik", logLik(fa), logLik(fb))
}
# cens event expression
fa <- tr("cens(ifelse(x > 1, 'right', 'none'))",
         frm(bf(y | cens(ifelse(x > 1, "right", "none")) ~ x), data = d))
d$cev <- ifelse(d$x > 1, "right", "none")
fb <- tr("cens(cev)", frm(bf(y | cens(cev) ~ x), data = d))
if (!is.null(fa)) same("cens event expr logLik", logLik(fa), logLik(fb))
# NA row in hi's variable on an interval row with na.omit, x NA
dcn <- dc; dcn$x[2] <- NA
fa <- tr("cens(cc, hi*1.5) with x NA", frm(bf(y_c | cens(cc, hi * 1.5) ~ x),
                                           data = dcn))
fb <- tr("cens(cc, hi2) with x NA", frm(bf(y_c | cens(cc, hi2) ~ x),
                                        data = dcn))
if (!is.null(fa)) same("cens interval expr, x NA logLik", logLik(fa),
                       logLik(fb))

cat("\n== 7. brms standata on each shape ==\n")
shapes <- list(
  list("weights(wt*2)", y | weights(wt * 2) ~ x, gaussian(),
       stats::gaussian(), "weights", "weights"),
  list("weights(wt,scale)", y | weights(wt, scale = TRUE) ~ x, gaussian(),
       stats::gaussian(), "weights", "weights"),
  list("se(s/2)", y | se(s / 2) ~ x, gaussian(), stats::gaussian(), "se",
       "se"),
  list("rate(t*2)", yc | rate(t * 2) ~ x, poisson(), stats::poisson(),
       "rate", "denom"),
  list("trials(c+1)", yb | trials(c + 1) ~ x, binomial(),
       stats::binomial(), "trials", "trials"),
  list("trials(max(c)+1)", yb | trials(max(c) + 1) ~ x, binomial(),
       stats::binomial(), "trials", "trials"),
  list("trunc(lb=min(y)-1)", y | trunc(lb = min(y) - 1) ~ x, gaussian(),
       stats::gaussian(), "trunc_lb", "lb"),
  list("trunc(ub=max(y)*2)", y | trunc(ub = max(y) * 2) ~ x, gaussian(),
       stats::gaussian(), "trunc_ub", "ub"),
  list("weights(log(wt)+1)", y | weights(log(wt) + 1) ~ x, gaussian(),
       stats::gaussian(), "weights", "weights"),
  list("weights(wt^2)", y | weights(wt^2) ~ x, gaussian(),
       stats::gaussian(), "weights", "weights"),
  list("se(sqrt(s))", y | se(sqrt(s)) ~ x, gaussian(), stats::gaussian(),
       "se", "se"),
  list("weights(1 + 0 * wt)", y | weights(1 + 0 * wt) ~ x, gaussian(),
       stats::gaussian(), "weights", "weights"))
for (sh in shapes) {
  fr <- tr(sh[[1]], frm(bf(sh[[2]]), data = d, family = sh[[3]],
                        dry_run = "frame"))
  s <- tr(paste("brms", sh[[1]]), sd_(sh[[2]], data = d, family = sh[[4]]))
  if (!is.null(fr) && !is.null(s)) {
    same(paste("standata", sh[[1]]), fr$aterm_values[[1]][[sh[[5]]]],
         as.numeric(s[[sh[[6]]]]))
  }
}
# a constant weight (all.vars empty): literal path
fr <- tr("weights(2)", frm(bf(y | weights(2) ~ x), data = d))
if (!is.null(fr)) cat("  weights(2) values head", head(fv(fr)$weights), "\n")
s <- tr("brms weights(2)", sd_(y | weights(2) ~ x, data = d))
if (!is.null(s)) cat("  brms weights(2) head", head(s$weights), "\n")
fr <- tr("weights(1 + 1)", frm(bf(y | weights(1 + 1) ~ x), data = d))
if (!is.null(fr)) cat("  weights(1+1) values head", head(fv(fr)$weights), "\n")

cat("\n== 8. readers of the frame on an expression model ==\n")
d$w2 <- d$wt * 2
fa <- frm(bf(y | weights(wt * 2) + trunc(lb = min(y) - 1) ~ x), data = d)
d$lbm <- min(d$y) - 1
fb <- frm(bf(y | weights(w2) + trunc(lb = lbm) ~ x), data = d)
same("readers logLik", logLik(fa), logLik(fb))
tr("residuals", same("residuals", residuals(fa), residuals(fb)))
tr("influence", same("influence cooks", influence(fa, groups = "f")$cooks,
                     influence(fb, groups = "f")$cooks))
tr("refit", same("refit logLik",
                 logLik(refit(fa, newresp = simulate(fa, seed = 1)[[1]])),
                 logLik(refit(fb, newresp = simulate(fb, seed = 1)[[1]]))))
tr("bootstrap", same("bootstrap t", frm_bootstrap(fa, nboot = 3, seed = 2)$t,
                     frm_bootstrap(fb, nboot = 3, seed = 2)$t))
tr("log_lik", same("log_lik", log_lik(fa), log_lik(fb)))
tr("loo", same("loo", loo(fa)$estimates, loo(fb)$estimates))
tr("simulate", same("simulate", simulate(fa, nsim = 2, seed = 3),
                    simulate(fb, nsim = 2, seed = 3)))
tr("pp_check data", {
  pa <- pp_check(fa, ndraws = 5); cat("  pp_check ok\n") })
tr("update(newdata)", same("update newdata logLik",
                           logLik(update(fa, newdata = d[1:40, ])),
                           logLik(update(fb, newdata = d[1:40, ]))))
cat("  names(fa$data):", names(fa$data), "\n")
cat("  names(fb$data):", names(fb$data), "\n")
# kfold refits read the frame
tr("kfold", same("kfold elpd", kfold(fa, K = 3, seed = 4)$estimates,
                 kfold(fb, K = 3, seed = 4)$estimates))
