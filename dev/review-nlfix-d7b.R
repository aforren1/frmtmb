# Fix round 2, item 7: does the early stop lose fits? REV_MUT = no_early restarts through all
# patterns; "" is the shipped rule. Reference: best of random starts. One family per process.
fam <- commandArgs(trailingOnly = TRUE)[1]
Sys.setenv(FRM_ARM = "mutant", REV_LIB = "C:/Users/adf44/source/r/rev-nlfix-lib/round3")
source("dev/review-nlfix-common.R")
q <- function(expr) suppressWarnings(suppressMessages(expr))
ll <- function(expr) { r <- tryCatch(q(expr), error = function(e) NULL); if (is.null(r)) NA_real_ else { v <- as.numeric(logLik(r)); if (is.finite(v)) v else NA_real_ } }
F <- list(cumulative = cumulative(), sratio = sratio(), cratio = cratio(), acat = acat(), hurdle_cumulative = hurdle_cumulative(), bernoulli = bernoulli(), poisson = poisson(), gaussian = gaussian(), gamma = Gamma("log"))
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
    cumulative = , sratio = , cratio = , acat = ordered(cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4)),
    hurdle_cumulative = ordered(ifelse(runif(n) < 0.2, 0L, as.integer(cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))), levels = 0:4),
    bernoulli = rbinom(n, 1, plogis(lat - 2)), poisson = rpois(n, lat), gaussian = lat + rnorm(n, 0, 0.5), gamma = rgamma(n, 3, 3 / lat))
  d
}
cat("family seed slope body default no_early ref t_def t_all\n")
for (seed in 1:5) for (slope in c(0.6, -0.6)) for (bn in names(bodies)) {
  b <- bodies[[bn]]; d <- gen(seed * 100 + 1, slope)
  fm <- eval(bquote(bf(y ~ .(b$f), nl = TRUE)))
  fm <- do.call(bf, c(list(as.formula(paste("y ~", deparse1(b$f)))), lapply(b$p, function(p) as.formula(paste(p, "~ 1"))), list(nl = TRUE)))
  Sys.setenv(REV_MUT = ""); t0 <- proc.time()[["elapsed"]]; a <- ll(frm(fm, family = F[[fam]], data = d)); t1 <- proc.time()[["elapsed"]]
  Sys.setenv(REV_MUT = "alias_only"); t2 <- proc.time()[["elapsed"]]; e <- ll(frm(fm, family = F[[fam]], data = d)); t3 <- proc.time()[["elapsed"]]
  Sys.setenv(REV_MUT = "")
  set.seed(seed * 7 + 1); best <- -Inf
  for (i in 1:8) {
    st <- setNames(rnorm(length(b$p), 0, 1), paste0(b$p, "_(Intercept)"))
    v <- ll(frm(fm, family = F[[fam]], data = d, start = list(beta = st)))
    if (!is.na(v) && v > best) best <- v
  }
  cat(sprintf("%s %d %.1f %s %.5f %.5f %.5f %.2f %.2f\n", fam, seed, slope, bn, a, e, best, t1 - t0, t3 - t2)); flush.console()
}
