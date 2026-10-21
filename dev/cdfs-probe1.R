# Which RTMB CDF primitives tape, and how good are their derivatives?
# Seed 1. Output: dev/cdfs-out/probe1.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
set.seed(1)
try_tape <- function(lab, f, x0) {
  r <- tryCatch({
    F <- MakeTape(f, x0)
    v <- F(x0)
    g <- F$jacobian(x0)
    h <- 1e-6 * pmax(1, abs(x0))
    fd <- vapply(seq_along(x0), function(j) {
      e <- replace(numeric(length(x0)), j, h[j])
      (f(x0 + e) - f(x0 - e)) / (2 * h[j])
    }, 0)
    sprintf("value %.15g  ref %.15g  grad %s  fd %s", v, f(x0),
            paste(signif(g, 8), collapse = ","),
            paste(signif(fd, 8), collapse = ","))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-40s %s\n", lab, r))
}
# pgamma: shape and scale on the tape
try_tape("pgamma lower", function(p) RTMB::pgamma(2.5, shape = exp(p[1]),
                                                   scale = exp(p[2])),
         c(log(3), log(0.7)))
try_tape("pgamma upper log", function(p) {
  RTMB::pgamma(2.5, shape = exp(p[1]), scale = exp(p[2]),
               lower.tail = FALSE, log.p = TRUE)
}, c(log(3), log(0.7)))
try_tape("pbeta", function(p) RTMB::pbeta(0.3, exp(p[1]), exp(p[2])),
         c(log(2), log(5)))
try_tape("pbeta upper", function(p) {
  RTMB::pbeta(0.3, exp(p[1]), exp(p[2]), lower.tail = FALSE)
}, c(log(2), log(5)))
try_tape("pbinom prob", function(p) RTMB::pbinom(3, 10, plogis(p[1])),
         c(-0.4))
try_tape("pnbinom mu", function(p) {
  RTMB::pnbinom(3, size = exp(p[1]), mu = exp(p[2]))
}, c(log(2), log(4)))
try_tape("pnbinom prob", function(p) {
  RTMB::pnbinom(3, size = exp(p[1]), prob = plogis(p[2]))
}, c(log(2), 0.3))
try_tape("ppois", function(p) RTMB::ppois(3, exp(p[1])), c(log(4)))
try_tape("pnorm log upper", function(p) {
  RTMB::pnorm(p[1], lower.tail = FALSE, log.p = TRUE)
}, c(1.3))
try_tape("pnorm log lower", function(p) {
  RTMB::pnorm(p[1], log.p = TRUE)
}, c(-1.3))
try_tape("sign", function(p) sign(p[1]) * p[1]^2, c(-0.7))
try_tape("besselI nu=3", function(p) besselI(exp(p[1]), 3), c(log(2.2)))
try_tape("besselI expon", function(p) besselI(exp(p[1]), 3, TRUE),
         c(log(2.2)))
try_tape("lgamma", function(p) lgamma(exp(p[1])), c(log(2.2)))
try_tape("log1p", function(p) log1p(p[1]), c(0.3))
try_tape("expm1", function(p) expm1(p[1]), c(0.3))
try_tape("RTMBdist::pt", function(p) RTMBdist::pt(p[1], exp(p[2])),
         c(0.8, log(5)))
try_tape("RTMBdist::pskewnorm", function(p) {
  RTMBdist::pskewnorm(p[1], 0, exp(p[2]), p[3])
}, c(0.4, 0.1, 2))
try_tape("RTMBdist::pexgauss", function(p) {
  RTMBdist::pexgauss(p[1], 0, exp(p[2]), exp(p[3]))
}, c(0.4, 0.1, 0.3))
try_tape("RTMBdist::pvm", function(p) {
  RTMBdist::pvm(1, p[1], exp(p[2]))
}, c(0.4, 0.3))
try_tape("pgamma upper (ppois tail)", function(p) {
  RTMB::pgamma(exp(p[1]), shape = 4)
}, c(log(4)))
