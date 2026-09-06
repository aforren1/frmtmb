library(RTMB)
tryfun <- function(nm, f) {
  r <- tryCatch({
    tp <- RTMB::MakeTape(function(p) sum(f(p)), c(0.7))
    v <- tp(0.7); g <- as.numeric(tp$jacobian(0.7))
    sprintf("OK  val=%.10g grad=%.10g", v, g)
  }, error = function(e) paste("ERR", conditionMessage(e)))
  cat(sprintf("%-16s %s\n", nm, r))
}
tryfun("pnorm", function(x) RTMB::pnorm(x))
tryfun("qnorm", function(x) RTMB::qnorm(x))
tryfun("atan", function(x) atan(x))
tryfun("asinh", function(x) asinh(x))
tryfun("sqrt", function(x) sqrt(x))
tryfun("expm1", function(x) expm1(x))
tryfun("log1p", function(x) log1p(x))
tryfun("abs", function(x) abs(x))
tryfun("logspace_add", function(x) RTMB::logspace_add(0*x, x))
tryfun("logspace_sub", function(x) RTMB::logspace_sub(x, 0*x))
tryfun("tan", function(x) tan(x))
tryfun("dnorm", function(x) dnorm(x))
cat("exists RTMB::qnorm:", exists("qnorm", envir = asNamespace("RTMB")), "\n")
cat("exists RTMB::pnorm:", exists("pnorm", envir = asNamespace("RTMB")), "\n")
# pnorm derivative against numDeriv
for (e in c(-3, -0.5, 0, 1.2, 4)) {
  tp <- RTMB::MakeTape(function(p) sum(RTMB::pnorm(p)), e)
  g <- as.numeric(tp$jacobian(e))
  nd <- numDeriv::grad(function(z) stats::pnorm(z), e)
  cat(sprintf("pnorm' at %5.2f  AD=%.16g  numDeriv=%.16g  dnorm=%.16g\n",
              e, g, nd, stats::dnorm(e)))
}
