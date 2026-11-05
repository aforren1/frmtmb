# Item 1: frm() against the closed form's global maximum, on designs the lane
# did not use. Args: lane|base  set  seeds.  Seeds in the design list below.
# Rscript dev/review-ncgp-truth.R lane fun 1:2
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- args[1]
set <- args[2]
seeds <- eval(parse(text = args[3]))
rev_libs(arm)
funs <- list(
  mix = function(x) 2 * sin(3 * x) + x,
  bump = function(x) 3 * exp(-(x - 5)^2),
  step = function(x) as.numeric(x > 5),
  line = function(x) 3 * x,
  cube = function(x) x^3 / 100,
  chirp = function(x) cos(x^2 / 5),
  logx = function(x) log1p(x)
)
row <- function(lab, y, x, seed) {
  d <- data.frame(y = y, x = x)
  t0 <- Sys.time()
  r <- rev_fit(y ~ gp(x), d)
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r$fit, "error")) {
    cat(sprintf("RT %s seed %d ERROR %s\n", lab, seed, conditionMessage(r$fit)))
    return(invisible())
  }
  g <- rev_gp_read(r$fit)
  s <- cf_setup(y, x, nug = if (arm == "lane") 1e-12 else 1e-6)
  xt <- list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd))))
  m <- cf_max(s, extra = xt)
  at <- cf_ll_at(s, g$beta, g$sigma, g$sd, g$rho)
  cat(sprintf(paste0("RT %s seed %d n %d m %d frm ll %.6f sigma %.5g sd %.5g rho ",
                     "%.5g conv %d nwarn %d | cf ll %.6f sigma %.5g sd %.5g rho ",
                     "%.5g | gap %.3e tape-vs-cf-at-point %.2e fit %.1fs\n"),
              lab, seed, length(y), s$m, g$ll, g$sigma, g$sd, g$rho, g$conv,
              length(r$warn), m$ll, m$sigma, m$sd, m$rho, g$ll - m$ll, g$ll - at, el))
  if (length(r$warn)) cat("RTW", substr(r$warn[1], 1, 110), "\n")
  flush.console()
}
if (set == "fun") {
  for (seed in seeds) for (fn in names(funs)) for (sg in c(1, 0.1, 1e-2, 1e-3)) {
    for (n in c(25L, 80L)) {
      set.seed(1000 * seed + n)
      x <- sort(round(runif(n, 0, 10), 4))
      y <- funs[[fn]](x) + sg * rnorm(n)
      row(sprintf("%s s=%g", fn, sg), y, x, seed)
    }
  }
}
if (set == "rep") {
  for (seed in seeds) for (mk in list(c(10, 5), c(20, 30), c(30, 10), c(15, 40), c(40, 3))) {
    for (sg in c(0.1, 0.02, 0.01, 3e-3, 1e-3)) {
      set.seed(7000 + seed)
      u <- sort(round(runif(mk[1], 0, 6), 3))
      x <- rep(u, each = mk[2])
      y <- sin(2 * x) + sg * rnorm(length(x))
      row(sprintf("rep %dx%d s=%g", mk[1], mk[2], sg), y, x, seed)
    }
  }
}
if (set == "dim2") {
  for (seed in seeds) for (sg in c(0.3, 0.03, 3e-3)) for (n in c(30L, 70L)) {
    set.seed(9000 + seed + n)
    x <- cbind(runif(n, 0, 5), runif(n, 0, 5))
    y <- sin(x[, 1]) * cos(x[, 2]) + sg * rnorm(n)
    d <- data.frame(y = y, x1 = x[, 1], x2 = x[, 2])
    r <- rev_fit(y ~ gp(x1, x2), d)
    if (inherits(r$fit, "error")) { cat("RT dim2 ERROR\n"); next }
    g <- rev_gp_read(r$fit)
    s <- cf_setup(y, x, nug = if (arm == "lane") 1e-12 else 1e-6)
    m <- cf_max(s, extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))))
    cat(sprintf("RT dim2 s=%g seed %d n %d frm ll %.6f sigma %.5g sd %.5g rho %.5g conv %d nwarn %d | cf ll %.6f sigma %.5g sd %.5g | gap %.3e\n",
                sg, seed, n, g$ll, g$sigma, g$sd, g$rho, g$conv, length(r$warn), m$ll,
                m$sigma, m$sd, g$ll - m$ll))
    if (length(r$warn)) cat("RTW", substr(r$warn[1], 1, 110), "\n")
  }
}
cat("DONE\n")
