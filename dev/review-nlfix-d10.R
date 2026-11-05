# Fix round 2, item 7: the shipped rule against the proposed one, on my data. REV_MUT:
# "" shipped, alias_only (no agreement stop; a first agreement stops only an exact alias),
# no_early, union. One data block per argument, min of 2.
blk <- commandArgs(trailingOnly = TRUE)[1]
Sys.setenv(FRM_ARM = "mutant", REV_LIB = "C:/Users/adf44/source/r/rev-nlfix-lib/round3")
source("dev/review-nlfix-common.R")
q <- function(expr) suppressWarnings(suppressMessages(expr))
arms <- c("", "alias_only", "no_early", "sign4")
go <- function(lab, f, k = 2) {
  for (m in arms) {
    Sys.setenv(REV_MUT = m); t <- Inf; r <- NULL
    for (i in seq_len(k)) { t <- min(t, system.time(r <- tryCatch(q(f()), error = function(e) e))[["elapsed"]]) }
    cat(sprintf("%-34s %-11s %7.2f s  ll %s code %s\n", lab, if (nzchar(m)) m else "shipped", t, if (inherits(r, "error")) conditionMessage(r) else format(as.numeric(logLik(r)), digits = 10), if (inherits(r, "error")) NA else r$opt$convergence)); flush.console()
  }
  Sys.setenv(REV_MUT = "")
}
if (blk == "cost") {
  set.seed(6401)
  mk <- function(n, ng) { d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(ng, n, TRUE))); d$y <- 2 * exp(0.3 * d$x) + rnorm(ng, 0, 0.3)[as.integer(d$g)] + rnorm(n, 0, 0.5); d$yo <- cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE); d }
  invisible(mk(2e5, 200)); invisible(mk(2e4, 200)); d3 <- mk(3e4, 30)
  fa <- bf(yo ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
  go("my c3 data, saddle, 3e4 rows", function() frm(fa, family = cumulative(), data = d3), k = 1)
  go("  an ordinary body, 3e4", function() frm(bf(yo ~ b1 * x, b1 ~ 1, nl = TRUE), family = cumulative(), data = d3), k = 2)
  set.seed(6403); n <- 2e4; ng <- 200
  d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(sample(ng, n, TRUE)))
  d$y <- 2 + 0.5 * d$x + rnorm(ng, 0, 0.3)[as.integer(d$g)] + rnorm(n, 0, 0.5)
  d$yo <- cut(0.8 * d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = 1:4, ordered_result = TRUE)
  go("RE, ridge a + b, 2e4", function() frm(bf(y ~ a + b + c0 * x, a ~ 1 + (1 | g), b ~ 1, c0 ~ 1, nl = TRUE), data = d))
  go("cumulative, ridge b0 + b1 x, 2e4", function() frm(bf(yo ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE), family = cumulative(), data = d))
}
if (blk == "lost") {
  gen <- function(seed, slope) { set.seed(seed); n <- 300; d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3)); lat <- 1.5 * exp(slope * d$x) + 0.3 * d$z; d$y <- ordered(cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4)); d }
  f1 <- bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE); f2 <- bf(y ~ b1 * exp(b2 * x) + b3 * z, b1 + b2 + b3 ~ 1, nl = TRUE)
  for (cs in list(list("acat", 101, f1), list("acat", 401, f1), list("acat", 501, f2), list("cumulative", 501, f2))) {
    d <- gen(cs[[2]], -0.6)
    for (m in arms) {
      Sys.setenv(REV_MUT = m)
      r <- tt(frm(cs[[3]], family = get(cs[[1]])(), data = d))
      cat(sprintf("%-11s seed %d %-11s ll %.4f code %s beta %s\n   M[%s]\n   W[%s]\n", cs[[1]], cs[[2]], if (nzchar(m)) m else "shipped", as.numeric(logLik(r)), r$opt$convergence, paste(format(r$estimates$beta, digits = 3), collapse = ","), substr(paste(attr(r, "m"), collapse = "|"), 1, 150), substr(paste(attr(r, "w"), collapse = "|"), 1, 150)))
    }
    Sys.setenv(REV_MUT = "")
  }
  cat("\nthe short rows no rule reaches: are they told?\n")
  d <- { set.seed(101); n <- 300; x <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3)); lat <- 1.5 * exp(0.6 * x$x) + 0.3 * x$z; x$y <- lat + rnorm(n, 0, 0.5); x }
  r <- tt(frm(bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE), data = d))
  cat(sprintf("gaussian sat: ll %.4f code %s beta %s\n M[%s]\n W[%s]\n", as.numeric(logLik(r)), r$opt$convergence, paste(format(r$estimates$beta, digits = 4), collapse = ","), substr(paste(attr(r, "m"), collapse = "|"), 1, 200), substr(paste(attr(r, "w"), collapse = "|"), 1, 260)))
  r2 <- tt(frm(bf(y ~ b1 * (1 - exp(-b2 * x)), b1 + b2 ~ 1, nl = TRUE), data = d, start = list(beta = c("b1_(Intercept)" = -2, "b2_(Intercept)" = -0.5))))
  cat(sprintf("  from b1 -2, b2 -0.5: ll %.4f code %s beta %s\n", as.numeric(logLik(r2)), r2$opt$convergence, paste(format(r2$estimates$beta, digits = 4), collapse = ",")))
  cat("\nmu body beside a disc body: union\n")
  set.seed(11001); n <- 300
  dd <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3)); invisible(factor(rep(1:15, each = 20))); invisible(rnorm(15, 0, 0.5))
  lat <- 1.5 * exp(0.6 * dd$x) + 0.7 * dd$z; dd$yo <- ordered(cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  fm4 <- bf(yo ~ 1) + nlf(mu ~ b1 * exp(b2 * x)) + nlf(disc ~ exp(d1 * d2 * z)) + lf(b1 ~ 1) + lf(b2 ~ 1) + lf(d1 ~ 1) + lf(d2 ~ 1)
  for (m in c("", "union", "sign4")) { Sys.setenv(REV_MUT = m); r <- tt(frm(fm4, family = cumulative(), data = dd)); cat(sprintf("%-8s ll %.4f code %s\n", if (nzchar(m)) m else "shipped", as.numeric(logLik(r)), r$opt$convergence)) }
  Sys.setenv(REV_MUT = "")
  r <- tt(frm(fm4, family = cumulative(), data = dd, start = list(beta = c("b1_(Intercept)" = 1, "b2_(Intercept)" = 0.5), betad = c(0.3, 0.3)))); cat(sprintf("all moved ll %.4f code %s\n", as.numeric(logLik(r)), r$opt$convergence))
}
