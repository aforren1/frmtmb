# The 450 nl rows behind "16 of 450 short" (dev/review-nlfix-d7.R): 9
# families, 5 bodies, slopes of both signs, 5 seeds. Each default fit
# against the best of 8 random starts, both on the arm's build. Records
# the default's warning classes and messages, so a short fit can be told
# silent from loud, and the wall time of the default fit.
# Copy of dev/nlstart-scan.R for lane nlstop (rows unchanged).
# Usage: NLSTOP_ARM=before|after [NLSTOP_REF=0] Rscript dev/nlstop-scan.R <family>
# Output: one line per row, prefixed ROW.
fam <- commandArgs(trailingOnly = TRUE)[1]
source("dev/nlstop-common.R")
llv <- function(r) {
  if (iserr(r) || is.null(r$opt)) return(NA_real_)
  v <- as.numeric(logLik(r))
  if (is.finite(v)) v else NA_real_
}
F <- list(cumulative = cumulative(), sratio = sratio(), cratio = cratio(),
          acat = acat(), hurdle_cumulative = hurdle_cumulative(),
          bernoulli = bernoulli(), poisson = poisson(), gaussian = gaussian(),
          gamma = Gamma("log"))
bodies <- list(
  exp = list(f = quote(b1 * exp(b2 * x)), p = c("b1", "b2")),
  sat = list(f = quote(b1 * (1 - exp(-b2 * x))), p = c("b1", "b2")),
  logi = list(f = quote(b1 / (1 + exp(b2 * x))), p = c("b1", "b2")),
  prod = list(f = quote(a + b1 * b2 * x), p = c("a", "b1", "b2")),
  two = list(f = quote(b1 * exp(b2 * x) + b3 * z), p = c("b1", "b2", "b3")))
gen <- function(seed, slope) {
  set.seed(seed); n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  lat <- 1.5 * exp(slope * d$x) + 0.3 * d$z
  d$y <- switch(fam,
    cumulative = , sratio = , cratio = , acat =
      ordered(cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4)),
    hurdle_cumulative = ordered(ifelse(runif(n) < 0.2, 0L, as.integer(
      cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))),
      levels = 0:4),
    bernoulli = rbinom(n, 1, plogis(lat - 2)), poisson = rpois(n, lat),
    gaussian = lat + rnorm(n, 0, 0.5), gamma = rgamma(n, 3, 3 / lat))
  d
}
cat("ROWHEAD family seed slope body default ref conv t_def nwarn wclass",
    "nmsg\n")
for (seed in 1:5) for (slope in c(0.6, -0.6)) for (bn in names(bodies)) {
  b <- bodies[[bn]]
  d <- gen(seed * 100 + 1, slope)
  fm <- do.call(bf, c(list(as.formula(paste("y ~", deparse1(b$f)))),
                      lapply(b$p, function(p) as.formula(paste(p, "~ 1"))),
                      list(nl = TRUE)))
  t0 <- proc.time()[["elapsed"]]
  r <- tt(frm(fm, family = F[[fam]], data = d))
  t1 <- proc.time()[["elapsed"]]
  a <- llv(r)
  conv <- if (iserr(r)) NA else r$opt$convergence
  wcl <- attr(r, "wc")
  set.seed(seed * 7 + 1)
  best <- -Inf
  for (i in seq_len(if (Sys.getenv("NLSTOP_REF", "1") == "0") 0 else 8)) {
    st <- setNames(rnorm(length(b$p), 0, 1), paste0(b$p, "_(Intercept)"))
    v <- llv(tt(frm(fm, family = F[[fam]], data = d,
                    start = list(beta = st))))
    if (!is.na(v) && v > best) best <- v
  }
  cat(sprintf("ROW %s %d %.1f %s %.6f %.6f %s %.3f %d %s %d\n", fam, seed,
              slope, bn, a, best, conv, t1 - t0, length(wcl),
              if (length(wcl)) paste(unique(wcl), collapse = "+") else "-",
              length(attr(r, "m"))))
  flush.console()
}
