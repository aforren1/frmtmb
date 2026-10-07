# Reviewer of lane setier: separation_check() on designs the lane's study
# did not cover. Per case: the frm() warnings (classified), code, and
# whether glm() / the LP reference call it separated.
#   Rscript dev/setier-rev-sep2.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-setier-lib",
           "C:/Users/adf44/source/r/rellib-r6"),
  base = "C:/Users/adf44/source/r/rellib-r6")
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
ns <- asNamespace("frmtmb")
cat("arm", arm, find.package("frmtmb"), "\n")
kind <- function(m) {
  if (grepl("separate the outcomes", m)) return("SEP")
  if (grepl("^Standard errors are not available", m)) return("SE")
  if (grepl("^Optimizer did not report", m)) return("CONV")
  if (grepl("gradient", m)) return("GRAD")
  if (grepl("^Boundary", m)) return("BOUNDARY")
  "OTHER"
}
run <- function(lab, expr) {
  w <- character(); m <- character()
  f <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
  }), error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("%-44s ERROR %s\n", lab, substr(conditionMessage(f), 1, 90)))
    return(invisible(NULL))
  }
  se <- suppressWarnings(suppressMessages(fixef(f)))[, "Est.Error"]
  k <- vapply(c(w, m), kind, "")
  cat(sprintf("%-44s code %d | %-16s | finite SE %d/%d | lost %s\n", lab,
              f$opt$convergence, if (length(k)) paste(k, collapse = "+")
              else "none", sum(is.finite(se)), length(se),
              paste(names(ns$sdr_of(f)$se_lost), collapse = ",")))
  sw <- w[grepl("separate the outcomes", w)]
  if (length(sw)) cat("    ", substr(sw[1], 1, 200), "\n")
  invisible(f)
}
cat("\n== a group with all-0 responses, no fixed-effect separation\n")
for (s in 1:5) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:15, each = 12)), x = rnorm(180))
  d$y <- rbinom(180, 1, plogis(0.3 * d$x + rnorm(15, 0, 1.5)[d$g]))
  d$y[d$g == "1"] <- 0L
  run(sprintf("glmm, group 1 all 0, seed %d", s),
      frm(y ~ x + (1 | g), family = bernoulli(), data = d))
}
cat("\n== a group with all-0 responses and a large group sd (sd 4)\n")
for (s in 1:5) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:15, each = 12)), x = rnorm(180))
  d$y <- rbinom(180, 1, plogis(0.3 * d$x + rnorm(15, 0, 4)[d$g]))
  run(sprintf("glmm sd 4, seed %d (zero groups %d, one groups %d)", s,
              sum(tapply(d$y, d$g, sum) == 0),
              sum(tapply(d$y, d$g, mean) == 1)),
      frm(y ~ x + (1 | g), family = bernoulli(), data = d))
}
cat("\n== response codings, complete separation y = x > 0\n")
set.seed(4)
d <- data.frame(x = rnorm(60), z = rnorm(60))
d$yi <- as.integer(d$x > 0)
d$yl <- d$x > 0
d$yf <- factor(ifelse(d$x > 0, "yes", "no"))
run("integer 0/1", frm(yi ~ x + z, family = bernoulli(), data = d))
run("logical", frm(yl ~ x + z, family = bernoulli(), data = d))
run("factor", frm(yf ~ x + z, family = bernoulli(), data = d))
cat("\n== binomial with trials, weights, offsets\n")
d$n <- sample(2:6, 60, TRUE)
d$k <- ifelse(d$x > 0, d$n, 0L)
run("trials, complete", frm(k | trials(n) ~ x + z, family = binomial(),
                            data = d))
d$k2 <- ifelse(d$x > 0, d$n, 0L); d$k2[1] <- 1L
run("trials, one row with both outcomes",
    frm(k2 | trials(n) ~ x + z, family = binomial(), data = d))
d$w <- ifelse(seq_len(60) <= 3, 0, 1)
d$yw <- d$yi; d$yw[1:3] <- 1L - d$yw[1:3]
run("weights 0 on the 3 rows that break it",
    frm(yw | weights(w) ~ x + z, family = bernoulli(), data = d))
run("same rows with weight 1 (not separated)",
    frm(yw ~ x + z, family = bernoulli(), data = d))
d$off <- rnorm(60)
run("offset, complete", frm(yi ~ x + z + offset(off), family = bernoulli(),
                            data = d))
run("beta_binomial trials, complete",
    frm(k | trials(n) ~ x + z, family = beta_binomial(), data = d))
cat("\n== other families with separated categories (not checked)\n")
set.seed(5)
d2 <- data.frame(x = rnorm(90))
d2$yc <- cut(d2$x, c(-Inf, -0.5, 0.5, Inf), labels = FALSE)
run("cumulative, perfectly ordered", frm(yc ~ x, family = cumulative(),
                                         data = d2))
d2$ycat <- factor(c("a", "b", "c")[d2$yc])
run("categorical, perfectly sorted", frm(ycat ~ x, family = categorical(),
                                         data = d2))
cat("\n== probit and cloglog links\n")
run("probit complete", frm(yi ~ x + z, family = bernoulli("probit"),
                           data = d))
run("cloglog complete", frm(yi ~ x + z, family = bernoulli("cloglog"),
                            data = d))
