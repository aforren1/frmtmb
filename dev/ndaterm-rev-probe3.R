# Reviewer of lane ndaterm: the other class-wise simulators (grouped
# mixture, hmm), thres(gr = ) in the CE grid, and draws paths.
# Usage: Rscript dev/ndaterm-rev-probe3.R <arm>   (mine | r9)
# Log: dev/ndaterm-rev-out/probe3-<arm>.txt (redirected by the caller)
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "mine"
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ndaterm/dev/stan-cache")
suppressMessages({library(frmtmb); library(frmtmb.latent)
  library(frmtmb.sample)})
cat("arm", arm, "frmtmb", find.package("frmtmb"), "sample",
    find.package("frmtmb.sample"), "\n")
tr <- function(label, expr) {
  w <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      structure(paste("ERROR:", conditionMessage(e)), class = "err")
    }),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning")
    })
  if (length(w)) cat("  [", label, "] WARN:", substr(unique(w), 1, 120),
                     "\n")
  if (inherits(r, "err")) cat("  [", label, "]", substr(r, 1, 300), "\n")
  invisible(r)
}
ok <- function(r) !inherits(r, "err") && !is.null(r)
show <- function(label, x) {
  cat(sprintf("%-52s %s\n", label,
              paste(format(x, digits = 5), collapse = "  ")))
}

cat("\n== 1. mixture(binomial, binomial, groups = ~id)\n")
set.seed(7)
ng <- 40; m <- 6
d <- data.frame(id = rep(seq_len(ng), each = m), x = rnorm(ng * m),
                N = sample(c(2, 50), ng * m, TRUE))
cl <- rep(rbinom(ng, 1, 0.5), each = m)
d$y <- rbinom(nrow(d), d$N, ifelse(cl == 1, 0.2, 0.8))
fg <- tr("group mixture fit", frm(bf(y | trials(N) ~ 1),
  family = mixture(binomial(), binomial(), groups = ~id), data = d))
if (ok(fg)) {
  s <- tr("simulate", as.matrix(simulate(fg, nsim = 100, seed = 1)))
  if (ok(s)) show("group mixture simulate: share of draws > N",
                  mean(s > d$N))
}

cat("\n== 2. hmm(K = 2, binomial()) with per-row trials()\n")
set.seed(7)
d2 <- data.frame(id = rep(1:20, each = 30), t = rep(1:30, 20),
                 N = sample(c(2, 50), 600, TRUE))
st <- unlist(lapply(1:20, function(i) {
  s <- numeric(30); s[1] <- 1
  for (j in 2:30) s[j] <- if (runif(1) < 0.9) s[j - 1] else 3 - s[j - 1]
  s
}))
d2$y <- rbinom(600, d2$N, c(0.2, 0.8)[st])
fh <- tr("hmm binomial fit", frm(bf(y | trials(N) ~ 1),
  family = hmm(K = 2, binomial(), time = t, group = id), data = d2))
if (ok(fh)) {
  s <- tr("hmm simulate", as.matrix(simulate(fh, nsim = 50, seed = 1)))
  if (ok(s)) show("hmm simulate: share of draws > N", mean(s > d2$N))
}

cat("\n== 3. cumulative thres(gr = g) in conditional_effects()\n")
set.seed(7)
d3 <- data.frame(x = rnorm(300), g = gl(3, 100))
d3$y <- cut(d3$x + rnorm(300) + c(-1, 0, 1)[d3$g], c(-Inf, -1, 0, 1, Inf),
            labels = FALSE)
f3 <- frm(bf(y | thres(gr = g) ~ x), family = cumulative(), data = d3)
r <- tr("CE default", conditional_effects(f3, "x"))
if (ok(r)) cat("thres(gr = g) CE default ok, grid g:",
               as.character(unique(r[[1]]$g)), "\n")
r <- tr("CE predict", conditional_effects(f3, "x", method = "predict",
                                          ndraws = 51))
if (ok(r)) cat("thres(gr = g) CE predict ok\n")
f3b <- tr("thres(gr = g) + (1 | g) fit", frm(
  bf(y | thres(gr = g) ~ x + (1 | g)), family = cumulative(), data = d3))
if (ok(f3b)) {
  r <- tr("CE (1|g) predict", conditional_effects(f3b, "x",
                                                  method = "predict",
                                                  ndraws = 51))
  if (ok(r)) cat("thres(gr = g) + (1 | g) CE predict ok\n")
  r <- tr("CE (1|g) epred", conditional_effects(f3b, "x"))
  if (ok(r)) cat("thres(gr = g) + (1 | g) CE epred ok\n")
}

cat("\n== 4. frmtmb.sample draws: binomial and a binomial mixture\n")
set.seed(7)
n <- 200
db <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE))
db$y <- rbinom(n, db$N, plogis(0.3 + 0.5 * db$x))
fb <- frm(bf(y | trials(N) ~ x), family = binomial(), data = db)
dr <- suppressWarnings(suppressMessages(frm_sample(
  fb, chains = 1, iter = 300, warmup = 150, refresh = 0, seed = 7)))
nd <- data.frame(x = c(0, 0), N = c(1000, 2), y = c(560, 1))
set.seed(1)
r <- tr("predictive_error(newdata)", predictive_error(dr, newdata = nd))
if (ok(r)) show("predictive_error(newdata) col means (y - yrep)",
                colMeans(r))
r <- tr("residuals(newdata)", residuals(dr, newdata = nd))
if (ok(r)) show("residuals(draws, newdata) Estimate", r[, "Estimate"])
r <- tr("predictive_interval(newdata)",
        predictive_interval(dr, newdata = nd))
if (ok(r)) show("predictive_interval(newdata)", r)
r <- tr("pp_check(draws, newdata)", pp_check(dr, newdata = nd,
                                             type = "stat", ndraws = 20))
if (ok(r)) cat("pp_check(draws, newdata) ok\n")
r <- tr("posterior_predict in-sample", posterior_predict(dr, ndraws = 20))
if (ok(r)) show("in-sample draws > N share", mean(sweep(r, 2, db$N, ">")))
r <- tr("posterior_predict re_formula NA nd",
        posterior_predict(dr, newdata = nd, re_formula = NA))
if (ok(r)) show("posterior_predict(newdata, re_formula = NA) means",
                colMeans(r))
r <- tr("CE draws method predict N=1000", conditional_effects(
  dr, "x", method = "posterior_predict",
  conditions = data.frame(N = 1000))[[1]])
if (ok(r)) show("CE(draws, predict, N = 1000) at x~0",
                r$estimate__[which.min(abs(r$x))])
dmv <- db
dmv$g2 <- 1 + dmv$x + rnorm(n)
fmv <- frm(mvbf(bf(y | trials(N) ~ x), bf(g2 ~ x)),
           family = list(binomial(), gaussian()), data = dmv)
dmvs <- suppressWarnings(suppressMessages(frm_sample(
  fmv, chains = 1, iter = 300, warmup = 150, refresh = 0, seed = 7)))
r <- tr("mv draws in-sample", posterior_predict(dmvs, resp = "y",
                                                ndraws = 10))
if (ok(r)) cat("mv draws in-sample ok\n")
r <- tr("mv draws resp g2 newdata without N",
        posterior_predict(dmvs, newdata = data.frame(x = 0), resp = "g2"))
if (ok(r)) show("mv draws resp = g2, no N", mean(r))
r <- tr("mv draws resp y newdata", posterior_predict(dmvs, newdata = nd,
                                                     resp = "y"))
if (ok(r)) show("mv draws resp = y, newdata means", colMeans(r))
cat("done\n")
