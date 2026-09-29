## Probe 2 on the REFERENCE build: construct replicates that LOSE the
## top ordinal category, then follow every refit path.
lib <- Sys.getenv("FRMTMB_PROBE_LIB",
                  "C:/Users/adf44/source/r/rellib-r3")
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("frmtmb", as.character(packageVersion("frmtmb")), "from", lib, "\n")
say <- function(...) cat(..., "\n", sep = "")

## n = 40, four categories, the top one rare (tau_3 = 2.6 on the logit
## scale with a slope of 0.5 puts about 6 percent of the rows there)
mk <- function(seed, n = 40, tau = c(-0.6, 0.5, 2.6), slope = 0.5) {
  set.seed(seed)
  x <- rnorm(n)
  eta <- slope * x
  cp <- cbind(plogis(tau[1] - eta), plogis(tau[2] - eta),
              plogis(tau[3] - eta), 1)
  u <- runif(n)
  y <- 1L + rowSums(u > cp[, 1:3, drop = FALSE])
  data.frame(x = x, y = as.integer(y))
}
dd <- mk(202)
say("data seed 202, n = 40, table(y) = ",
    paste(table(factor(dd$y, levels = 1:4)), collapse = "/"))

fam_thresholds <- function(f) length(f$estimates$tau_raw)

for (famnm in c("cumulative", "sratio")) {
  f <- get(famnm, envir = asNamespace("frmtmb"))()
  fit <- suppressWarnings(frm(bf(y ~ x), family = f, data = dd))
  say("== ", famnm, ": fitted thresholds = ", fam_thresholds(fit),
      " logLik = ", format(as.numeric(logLik(fit)), digits = 12))

  set.seed(11)
  sims <- simulate(fit, nsim = 60, re_formula = NA)
  tops <- vapply(sims, function(v) max(as.integer(v)), 1L)
  lost <- which(tops < 4L)
  say(famnm, ": seed 11, nsim 60: replicates with max(y) < 4: ",
      length(lost), " (", paste(utils::head(lost, 12), collapse = ","), ")")
  if (!length(lost)) next

  b <- lost[1L]
  yb <- as.integer(sims[[b]])
  say(famnm, ": replicate ", b, " table = ",
      paste(table(factor(yb, levels = 1:4)), collapse = "/"))

  rf <- tryCatch(suppressWarnings(refit(fit, yb)),
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(rf)) {
    say(famnm, ": refit() -> ", rf)
  } else {
    say(famnm, ": refit() thresholds = ", fam_thresholds(rf),
        " logLik = ", format(as.numeric(logLik(rf)), digits = 15))
  }
  db <- data.frame(x = dd$x, y = yb)
  fp <- tryCatch(suppressWarnings(
    frm(bf(y | thres(3) ~ x), family = f, data = db)),
    error = function(e) paste("ERROR:", conditionMessage(e)))
  say(famnm, ": thres(3) pinned -> ",
      if (is.character(fp)) fp else
        paste0("thresholds = ", fam_thresholds(fp), " logLik = ",
               format(as.numeric(logLik(fp)), digits = 15)))
  fdd <- tryCatch(suppressWarnings(frm(bf(y ~ x), family = f, data = db)),
                  error = function(e) paste("ERROR:", conditionMessage(e)))
  say(famnm, ": default fresh -> ",
      if (is.character(fdd)) fdd else
        paste0("thresholds = ", fam_thresholds(fdd), " logLik = ",
               format(as.numeric(logLik(fdd)), digits = 15)))

  ## bootstrap with a FUN that reports the thresholds too
  FUNt <- function(ff) {
    c(fixef(ff, flatten = TRUE), tau = ff$estimates$tau_raw)
  }
  bs <- tryCatch(suppressWarnings(
    frm_bootstrap(fit, FUN = FUNt, nsim = 20, seed = 11)),
    error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(bs)) {
    say(famnm, ": frm_bootstrap(FUN incl thresholds) -> ", bs)
  } else {
    say(famnm, ": boot t dim = ", paste(dim(bs$t), collapse = "x"),
        " colnames = ", paste(colnames(bs$t), collapse = ","))
    say(famnm, ": all-NA rows = ",
        sum(apply(bs$t, 1, function(r) all(is.na(r)))),
        " any-NA rows = ", sum(apply(bs$t, 1, anyNA)))
  }
}

## ---- grouped thresholds -------------------------------------------
mkg <- function(seed, n = 90) {
  set.seed(seed)
  g <- factor(rep(c("a", "b", "c"), length.out = n))
  x <- rnorm(n)
  tau <- list(a = c(-0.6, 0.6, 2.4), b = c(-0.5, 1.2), c = c(-0.3, 0.9, 2.7))
  y <- integer(n)
  for (i in seq_len(n)) {
    tt <- tau[[as.character(g[i])]]
    cp <- c(plogis(tt - 0.5 * x[i]), 1)
    y[i] <- 1L + sum(runif(1) > cp[-length(cp)])
  }
  data.frame(x = x, g = g, y = y)
}
dg <- mkg(303)
say("grouped: seed 303 n = 90, per-level tables:")
print(table(dg$g, dg$y))
fg <- suppressWarnings(frm(bf(y | thres(gr = g) ~ x),
                          family = cumulative(), data = dg))
say("grouped: tau_raw length = ", length(fg$estimates$tau_raw))
set.seed(13)
sg <- simulate(fg, nsim = 60, re_formula = NA)
lostg <- vapply(sg, function(v) {
  yy <- as.integer(v)
  any(tapply(yy, dg$g, max) < tapply(dg$y, dg$g, max))
}, TRUE)
say("grouped: seed 13, nsim 60: replicates where some level lost its ",
    "top category: ", sum(lostg), " (",
    paste(utils::head(which(lostg), 12), collapse = ","), ")")
if (any(lostg)) {
  b <- which(lostg)[1L]
  yb <- as.integer(sg[[b]])
  print(table(dg$g, yb))
  rf <- tryCatch(suppressWarnings(refit(fg, yb)),
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  say("grouped: refit() -> ",
      if (is.character(rf)) rf else
        paste0("tau_raw length = ", length(rf$estimates$tau_raw),
               " logLik = ", format(as.numeric(logLik(rf)), digits = 15)))
  dbg <- data.frame(x = dg$x, g = dg$g, y = yb)
  npin <- tapply(dg$y, dg$g, max) - 1L
  say("grouped: fitted per-level counts = ",
      paste(npin, collapse = ","))
  ## pin per level through thres(n, gr = g)
  dbg$k <- as.integer(npin[as.character(dbg$g)])
  fpg <- tryCatch(suppressWarnings(
    frm(bf(y | thres(k, gr = g) ~ x), family = cumulative(), data = dbg)),
    error = function(e) paste("ERROR:", conditionMessage(e)))
  say("grouped: thres(k, gr = g) pinned -> ",
      if (is.character(fpg)) fpg else
        paste0("tau_raw length = ", length(fpg$estimates$tau_raw),
               " logLik = ", format(as.numeric(logLik(fpg)), digits = 15)))
  fdg <- tryCatch(suppressWarnings(
    frm(bf(y | thres(gr = g) ~ x), family = cumulative(), data = dbg)),
    error = function(e) paste("ERROR:", conditionMessage(e)))
  say("grouped: default fresh -> ",
      if (is.character(fdg)) fdg else
        paste0("tau_raw length = ", length(fdg$estimates$tau_raw),
               " logLik = ", format(as.numeric(logLik(fdg)), digits = 15)))
}

## ---- B, with a valid dummy response -------------------------------
ddb <- data.frame(x = rnorm(40), y = 1L)
for (famnm in c("cratio", "acat", "sratio", "cumulative")) {
  f <- get(famnm, envir = asNamespace("frmtmb"))()
  r <- tryCatch({
    s <- frm_simulate(bf(y | thres(3) ~ x), family = f, data = ddb,
                      prior = set_prior("normal(0, 2)",
                                        class = "Intercept") +
                        set_prior("normal(0, 1)", class = "b"),
                      nsim = 3, seed = 5)
    list(ok = TRUE, pars = attr(s, "pars"), y = s)
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  if (!r$ok) {
    say("B: ", famnm, " -> ERROR: ", r$msg)
  } else {
    say("B: ", famnm, " -> pars columns: ",
        paste(names(r$pars), collapse = ", "))
    print(r$pars)
    say("B: ", famnm, " simulated category counts per sim:")
    print(vapply(r$y, function(v) table(factor(as.integer(v),
                                               levels = 1:4)), integer(4)))
  }
}
