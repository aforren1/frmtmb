# Reviewer of lane setier, re-check: separated fits that sep_certificate()
# can miss after punch round 1 (the estimate's own direction first; a
# null space only when some tail is below 1e-3 and p <= 2000). Against
# the LP reference of dev/setier-sep.R.
#   Rscript dev/setier-rev2-sepmiss.R <lib> [part ...]
# Parts:
#   budget  complete separation (dev/setier-sep.R "complete", n = 50) and
#           factor quasi-separation ("quasi") under eval.max 5, 10, 20, 40
#           and the default, seeds 1..20
#   few     quasi-separation by a level with 2 rows, both 0, n = 40,
#           seeds 1..20
#   bigp    quasi-separation in a 2050-level factor (sparse_x), seed 1
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
parts <- if (length(args) > 1) args[-1] else c("budget", "few", "bigp")
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", find.package("frmtmb"), "\n")
src <- readLines("dev/setier-sep.R")
eval(parse(text = src[grep("^lp_separated <- function", src):
                        (grep("^gen <- function", src) - 1L)]))
kind <- function(w) {
  if (!length(w)) return("none")
  paste(vapply(w, function(m) {
    if (grepl("separate the outcomes", m)) "SEP" else
      if (grepl("^Optimizer did not report", m)) "CONV" else
        if (grepl("^Standard errors are not available", m)) "SE" else
          if (grepl("gradient", m)) "GRAD" else "OTHER"
  }, ""), collapse = "+")
}
fitw <- function(expr) {
  w <- character()
  f <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) invokeRestart("muffleMessage"))
  list(f = f, w = w)
}
tails <- function(f) {
  lp <- f$frame$linpreds[[1]]
  mu <- lp$link$linkinv(as.numeric(lp$X %*% f$estimates$beta[lp$idx]))
  y <- f$frame$y[[1]]
  min(ifelse(y > 0, 1 - mu, mu))
}
if ("budget" %in% parts) {
  cat("\n== budgets\n")
  for (des in c("complete", "quasi")) for (ev in c(5, 10, 20, 40, NA)) {
    named <- 0; sep_lp <- 0; codes <- 0; kinds <- character()
    for (s in 1:20) {
      set.seed(s)
      if (des == "complete") {
        d <- data.frame(x = rnorm(50), z = rnorm(50))
        d$y <- as.integer(d$x > 0); fo <- y ~ x + z
      } else {
        d <- data.frame(f = factor(rep(c("a", "b", "c"), each = 40)),
                        x = rnorm(120))
        d$y <- rbinom(120, 1, plogis(0.3 + 0.5 * d$x))
        d$y[d$f == "c"] <- 0L; fo <- y ~ f + x
      }
      sep_lp <- sep_lp + lp_separated(model.matrix(fo, d), d$y, 1)
      ctl <- if (is.na(ev)) frmtmb_control() else
        frmtmb_control(optCtrl = list(eval.max = ev, iter.max = ev),
                       restarts = 0)
      r <- fitw(frm(fo, family = bernoulli(), data = d, control = ctl))
      named <- named + any(grepl("separate the outcomes", r$w))
      codes <- codes + (r$f$opt$convergence != 0)
      kinds <- c(kinds, kind(r$w))
    }
    tb <- table(kinds)
    cat(sprintf("%-8s eval.max %-4s | LP separated %d | named %d | code!=0 %d | %s\n",
                des, if (is.na(ev)) "dflt" else ev, sep_lp, named, codes,
                paste(names(tb), tb, collapse = ", ")))
  }
}
if ("few" %in% parts) {
  cat("\n== a level with 2 rows, both 0\n")
  named <- 0; sep_lp <- 0; kinds <- character(); mint <- numeric()
  for (s in 1:20) {
    set.seed(s)
    d <- data.frame(f = factor(c(rep("a", 19), rep("b", 19), "c", "c")),
                    x = rnorm(40))
    d$y <- rbinom(40, 1, plogis(0.3 + 0.5 * d$x))
    d$y[d$f == "c"] <- 0L
    sep_lp <- sep_lp + lp_separated(model.matrix(y ~ f + x, d), d$y, 1)
    r <- fitw(frm(y ~ f + x, family = bernoulli(), data = d))
    named <- named + any(grepl("separate the outcomes", r$w))
    kinds <- c(kinds, kind(r$w))
    mint <- c(mint, tails(r$f))
  }
  tb <- table(kinds)
  cat(sprintf("LP separated %d | named %d | smallest tail median %.3g max %.3g | %s\n",
              sep_lp, named, median(mint), max(mint),
              paste(names(tb), tb, collapse = ", ")))
}
if ("bigp" %in% parts) {
  cat("\n== 2050-level factor, sparse_x\n")
  set.seed(1)
  L <- 2050
  d <- data.frame(f = factor(rep(seq_len(L), each = 4)), x = rnorm(4 * L))
  d$y <- rbinom(4 * L, 1, plogis(0.5 + 0.5 * d$x))
  zero <- tapply(d$y, d$f, sum) == 0
  one <- tapply(d$y, d$f, mean) == 1
  t0 <- proc.time()[["elapsed"]]
  r <- fitw(frm(y ~ f + x, family = bernoulli(), data = d,
                control = frmtmb_control(sparse_x = TRUE)))
  t1 <- proc.time()[["elapsed"]]
  cat(sprintf("levels all 0: %d, all 1: %d | fit %.1f s code %d | warnings: %s\n",
              sum(zero), sum(one), t1 - t0, r$f$opt$convergence, kind(r$w)))
  cat("  smallest tail:", signif(tails(r$f), 3), "\n")
  lost <- ns$sdr_of(r$f)$se_lost
  cat("  lost:", length(lost), "reasons:",
      paste(names(table(lost)), table(lost), collapse = ", "), "\n")
  cat("  first warning:", substr(c(r$w, "")[1], 1, 200), "\n")
}
