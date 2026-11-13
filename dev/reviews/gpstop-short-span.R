# Short check of lane gpstop's fix round: the length-scale span rule and
# the 1e-2 inner eigenvalue cut of gp_relin_regular(). Each call of the
# check prints its inputs (REG): code, gradient, length scale against the
# smallest gap, the MEDIAN nearest-neighbour gap and the range, the inner
# eigenvalue, the verdict. Arms lane and norelin interleaved.
# Sets (arg 1):
#   tie  - bernoulli bodies, n 40, seeds 1:6, positions as drawn ("plain")
#          and with one near-tie (x[20] <- x[19] + 1e-4, "tie")
#   long - gaussian bodies whose field is a linear trend in x
# Data: set.seed(82200 + 100 * seed + n + 10 * body index (+ 5 for long)).
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/reviews/gpstop-review-lib.R")
rv_libs("lane")
rv_instrument()
ns <- asNamespace("frmtmb")
real_reg <- get("gp_relin_regular", ns)
CNT$reg_real <- real_reg
rv_put("gp_relin_regular", function(obj, frame, o, tol) {
  v <- .GlobalEnv$CNT$reg_real(obj, frame, o, tol)
  pl <- obj$env$parList(o$par)
  th <- as.numeric(pl[["theta"]])
  out <- character(0)
  for (bk in frame[["re_blocks"]] %||% list()) {
    if (!identical(bk[["covstruct"]], "gp")) next
    d2 <- Reduce(`+`, bk[["aux_D2"]])
    d <- sqrt(d2[d2 > 0])
    D <- sqrt(d2)
    diag(D) <- Inf
    nn <- stats::median(apply(D, 1, min))
    ls <- exp(th[bk[["theta_idx"]][-1L]])
    out <- c(out, sprintf("ls %.3g ls/mingap %.3g ls/mednn %.3g ls/range %.3g",
                          ls, ls / min(d), ls / nn, ls / max(d)))
  }
  g <- max(abs(obj$gr(o$par)))
  ev <- tryCatch(gp_inner_ev_min(obj, o$par), error = function(e) NA)
  cat(sprintf("REG %s conv %d gmax %.2e %s inner_ev %.3e verdict %s\n",
              .GlobalEnv$CNT$tag, as.integer(o$convergence), g,
              paste(out, collapse = " | "), ev, v))
  v
})
bodies <- c("exp_a", "logis_a", "pow_a", "mm_V")
spell <- function(body) switch(body,
  exp_a = bf(y ~ a * exp(b * z), a ~ gp(x), b ~ 1, nl = TRUE),
  logis_a = bf(y ~ a / (1 + exp(-b * z)), a ~ gp(x), b ~ 1, nl = TRUE),
  pow_a = bf(y ~ a * w^b, a ~ gp(x), b ~ 1, nl = TRUE),
  mm_V = bf(y ~ V * w / (K + w), V ~ gp(x), K ~ 1, nl = TRUE))
mk <- function(body, n, seed, long = FALSE) {
  set.seed(82200 + 100 * seed + n + 10 * match(body, bodies) +
             if (long) 5 else 0)
  x <- sort(round(runif(n, 0, 10), 3))
  z <- round(rnorm(n), 3)
  w <- round(runif(n, 0.2, 5), 3)
  f <- if (long) 0.15 * x - 0.75 else 2 * sin(3 * x) + x
  fs <- if (long) f else f / max(abs(f))
  eta <- switch(body,
    exp_a = (2 + fs) * exp(0.4 * z),
    logis_a = (3 + fs) / (1 + exp(-1.5 * z)),
    pow_a = (2 + fs) * w^0.5,
    mm_V = (3 + fs) * w / (1 + w))
  y <- if (long) eta + 0.1 * rnorm(n) else
    rbinom(n, 1, plogis(1.5 * (eta - median(eta)) / sd(eta)))
  data.frame(x, z, w, y)
}
set <- a[1]
i <- 0L
if (set == "tie") {
  for (seed in 1:6) for (body in bodies) for (v in c("plain", "tie")) {
    d <- mk(body, 40L, seed)
    if (v == "tie") d$x[20] <- d$x[19] + 1e-4
    tag <- sprintf("%s_bern_n40_s%d_%s", body, seed, v)
    CNT$tag <- tag
    i <- i + 1L
    arms <- if (i %% 2) c("lane", "norelin") else c("norelin", "lane")
    for (arm in arms) rv_run(tag, arm, spell(body), data = d,
                             family = bernoulli())
  }
} else {
  for (seed in 1:3) for (n in c(60L, 120L)) for (body in c("exp_a", "pow_a")) {
    d <- mk(body, n, seed, long = TRUE)
    tag <- sprintf("%s_long_n%d_s%d", body, n, seed)
    CNT$tag <- tag
    i <- i + 1L
    arms <- if (i %% 2) c("lane", "norelin") else c("norelin", "lane")
    for (arm in arms) rv_run(tag, arm, spell(body), data = d)
  }
}
cat("DONE\n")
