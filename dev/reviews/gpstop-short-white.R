# Short check of lane gpstop's fix round: does gp_relin_regular() accept a
# white-noise field when one pair of positions is close? MAJOR 2's data
# (gpstop-review-sweep.R mk("logis_a", "bernoulli", 40, 2)) with x[20] moved
# to x[19] + 1e-4. The first build's degenerate point (beta -11.046,
# -135.46, theta 3.5988, -323.27) is placed on the tape with the log length
# scale set to log(k * mingap / 10), k = 2, 5, 20, and polished from there
# with nlminb; the rule is asked of the result.
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/reviews/gpstop-review-lib.R")
rv_libs("lane")
ns <- asNamespace("frmtmb")
bodies <- c("exp_a", "exp_b", "logis_a", "mm_V", "mm_K", "pow_a", "pow_b")
set.seed(81100 + 1000 * 3 + 100 * 2 + 40 + 10 * match("logis_a", bodies))
n <- 40L
x <- sort(round(runif(n, 0, 10), 3))
z <- round(rnorm(n), 3)
w <- round(runif(n, 0.2, 5), 3)
f <- 2 * sin(3 * x) + x
fs <- f / max(abs(f))
eta <- (3 + fs) / (1 + exp(-1.5 * z))
y <- rbinom(n, 1, plogis(1.5 * (eta - median(eta)) / sd(eta)))
d <- data.frame(x, z, w, y)
fml <- bf(y ~ a / (1 + exp(-b * z)), a ~ gp(x), b ~ 1, nl = TRUE)
for (v in c("plain", "tie")) {
  dd <- d
  if (v == "tie") dd$x[20] <- dd$x[19] + 1e-4
  q <- rv_fit(fml, data = dd, family = bernoulli())
  fit <- q$fit
  cat(sprintf("WHITE %s default ll %.6f conv %d nwarn %d theta %s\n", v,
              as.numeric(logLik(fit)), fit$opt$convergence, length(q$warn),
              paste(sprintf("%.3f", fit$opt$par[names(fit$opt$par) == "theta"]),
                    collapse = ",")))
  obj <- fit$obj
  frame <- fit$frame
  bk <- Filter(function(b) identical(b$covstruct, "gp"), frame$re_blocks)[[1]]
  d2 <- Reduce(`+`, bk$aux_D2)
  mg <- min(sqrt(d2[d2 > 0]))
  pn <- names(obj$par)
  for (k in c(2, 5, 20)) {
    p <- obj$par
    p[pn == "beta"] <- c(-11.046, -135.46)
    p[pn == "theta"] <- c(3.5988, log(k * mg / 10))
    o <- tryCatch(stats::nlminb(p, obj$fn, obj$gr,
                                control = list(iter.max = 200)),
                  error = function(e) NULL)
    if (is.null(o)) { cat("WHITE", v, k, "nlminb error\n"); next }
    o$convergence <- as.integer(o$convergence)
    reg <- tryCatch(ns$gp_relin_regular(obj, frame, o, 1e-3),
                    error = function(e) paste("ERR", conditionMessage(e)))
    g <- max(abs(obj$gr(o$par)))
    ev <- ns$gp_inner_ev_min(obj, o$par)
    cat(sprintf(paste0("WHITE %s k %g mingap %.3g start ls %.3g | end ll %.6f ",
                       "conv %d gmax %.2e inner_ev %.3g theta %s | ",
                       "ls/mingap %.3g | regular %s\n"),
                v, k, mg, k * mg / 10, -o$objective, o$convergence, g, ev,
                paste(sprintf("%.3f", o$par[pn == "theta"]), collapse = ","),
                exp(o$par[pn == "theta"][2]) / mg, format(reg)))
    flush.console()
  }
}
cat("DONE\n")
