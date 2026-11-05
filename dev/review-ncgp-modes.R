# Item 1: the silent short stops. What the start scan reads, and which cheap
# remedy reaches the marginal's maximum. Seeds below. Lane build, in process.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
# a coarse closed-form grid: rho from the largest distance down, r = sigma^2/sd^2
# from 1e-2 to 1e-10. The cost is 6 x 5 Cholesky factors of the position matrix.
cf_coarse <- function(y, x) {
  s <- cf_setup(y, x)
  d <- sqrt(max(s$D2))
  best <- NULL
  for (rho in d / 4^(0:5)) for (lr in log(10^c(1, -1, -3, -6, -9))) {
    v <- cf_prof(s, log(rho), lr)
    if (is.null(best) || (is.finite(v$ll) && v$ll > best$ll)) best <- v
  }
  best
}
arms <- list(
  lane = function(y, x) rev_fit(y ~ gp(x), data.frame(y = y, x = x)),
  rho1 = function(y, x) rev_fit(y ~ gp(x), data.frame(y = y, x = x), start = list(theta = c(0, 0))),
  cfstart = function(y, x) {
    t0 <- Sys.time()
    b <- cf_coarse(y, x)
    r <- rev_fit(y ~ gp(x), data.frame(y = y, x = x),
                 start = list(theta = c(log(b$sd), log(b$rho)), betad = log(b$sigma), beta = b$beta[1]))
    r
  }
)
rows <- list()
for (fn in c("mix", "step", "logx", "cube", "sin6", "sin6k")) for (sg in c(0.1, 1e-2, 1e-3)) for (n in c(25L, 80L)) {
  set.seed(1000 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  e <- rnorm(n)
  y <- switch(fn, mix = 2 * sin(3 * x) + x, step = as.numeric(x > 5), logx = log1p(x), cube = x^3 / 100,
              sin6 = sin(6 * x), sin6k = 1000 * sin(6 * x)) + sg * e * (if (fn == "sin6k") 1000 else 1)
  s <- cf_setup(y, x)
  m <- cf_max(s, ngrid = 30L)
  out <- sprintf("RM %-6s s=%g n %d cf ll %.4f sigma %.4g sd %.4g rho %.4g |", fn, sg, n, m$ll, m$sigma, m$sd, m$rho)
  for (a in names(arms)) {
    t0 <- Sys.time()
    r <- arms[[a]](y, x)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(r$fit, "error")) { out <- paste(out, a, "ERR |"); next }
    g <- rev_gp_read(r$fit)
    out <- paste0(out, sprintf(" %s gap %.2e conv %d w %d %.1fs |", a, g$ll - m$ll, g$conv, length(r$warn), el))
  }
  cat(out, "\n"); flush.console()
}
cat("== what the scan sees on mix\n")
set.seed(1080); x <- sort(round(runif(80, 0, 10), 4)); y <- 2 * sin(3 * x) + x + 1e-3 * rnorm(80)
f <- suppressWarnings(frm(y ~ gp(x), data.frame(y, x), control = NULL))
bk <- f$frame$re_blocks[[1]]
print(ns$gp_start(bk))
d <- sqrt(max(Reduce(`+`, bk$aux_D2)))
tp <- f$frame$par_template
cat("par_template betad, beta:", tp$betad, tp$beta, "\n")
ob <- f$obj
for (r in d / 4^(0:3)) {
  p <- ns$outer_from_template({t2 <- tp; t2$theta[2] <- log(r); t2}, ob, f$frame)
  cat(sprintf("scan rho %.4g nll %.4f\n", r, ob$fn(p)))
}
cat("DONE\n")
