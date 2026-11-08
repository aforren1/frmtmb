# Lane latperf: how often a fitted row's pair of positions leaves the pattern
# of the inner factor, so that jcs_rowvar_inner() solves it directly. Counts
# the rows each call of jcs_rowvar_solve() takes, on the designs of
# dev/latperf-agree.R (seed 11) and dev/latperf-prof.R (seed 1, 8000 rows).
.libPaths(c("C:/Users/adf44/source/r/wt-latperf-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
solved <- 0L
trace("jcs_rowvar_solve", tracer = quote(solved <<- solved + nrow(Ar)),
      where = ns, print = FALSE)
selinv <- 0L
trace("jcs_selinv_keys", tracer = quote(selinv <<- selinv + 1L),
      where = ns, print = FALSE)
src <- readLines("dev/latperf-agree.R")
# the generator and the cases of the agreement script, without its loop
eval(parse(text = src[grep("^mk <- ", src):(grep("^dense_copy", src) - 1L)]))
eval(parse(text = src[grep("^cases <- ", src):(grep("^d <- mk", src))]))
for (cs in cases) {
  fit <- suppressWarnings(suppressMessages(
    frm(stats::as.formula(cs[[2]]), family = cs[[3]], data = d)))
  solved <- 0L
  selinv <- 0L
  invisible(suppressWarnings(frm_linpred(fit, se.fit = TRUE)))
  cat(sprintf("RESULT %-22s rows %d  selected-inverse calls %d  rows solved %d\n",
              cs[[1]], nrow(d), selinv, solved))
}
for (design in c("one", "two")) {
  n <- 8000L
  set.seed(1)
  d2 <- data.frame(t = seq_len(n), g = factor(1), x = rnorm(n),
                   id = factor(sample.int(50, n, TRUE)))
  e <- as.numeric(arima.sim(list(ar = .5), n, sd = .5))
  u <- rnorm(50, 0, 0.3)
  d2$y <- rpois(n, exp(.5 + .3 * d2$x + e +
                        if (design == "two") u[d2$id] else 0))
  f <- if (design == "two") y ~ x + (1 | id) + ar(t, g, cov = TRUE) else
    y ~ x + ar(t, g, cov = TRUE)
  fit <- suppressWarnings(frm(f, family = poisson(), data = d2))
  solved <- 0L
  selinv <- 0L
  invisible(frm_linpred(fit, se.fit = TRUE))
  cat(sprintf("RESULT 8000 %-15s rows %d  selected-inverse calls %d  rows solved %d\n",
              design, n, selinv, solved))
}
