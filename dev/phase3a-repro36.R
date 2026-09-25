# Item 3.6: the all-censored group that rp_floored() cannot see.
#
# The construction is dev/frailty/frailty-floor3.R's, seed 20260910:
# 400 subjects, 40 centres of 10, population shape 0.6, centre sd 0.35,
# five centres followed to a common administrative time with no deaths.
# Run with PHASE3A_ARM=base for the released build and =lane for this
# lane's build. Seeds 20260910 to 20260915, as in the frailty lane.
source("C:/Users/adf44/source/r/frmtmb-wt-phase3a/dev/phase3a-lib.R")
suppressMessages({
  library(frmtmb); library(frmtmb.spline)
})
phase3a_where("frmtmb"); phase3a_where("frmtmb.spline")

sim_nodeath <- function(seed, n = 400L, n_centre = 40L, n_nodeath = 5L,
                        sd_u = 0.35, gamma1 = 0.6, scale = 5,
                        beta = 0.6, p_cens = 0.4, floor_true = 0.05) {
  set.seed(seed)
  centre <- rep(seq_len(n_centre), length.out = n)
  repeat {
    u <- stats::rnorm(n_centre, 0, sd_u)
    if (all(gamma1 + u > floor_true)) break
  }
  trt <- stats::rbinom(n, 1L, 0.5)
  g0 <- -gamma1 * log(scale)
  tt <- exp((log(-log(stats::runif(n))) - g0 - beta * trt) /
              (gamma1 + u[centre]))
  tau <- unname(stats::quantile(tt, 1 - p_cens))
  ev <- as.integer(tt <= tau)
  tt <- pmin(tt, tau)
  if (n_nodeath > 0L) {
    j <- centre <= n_nodeath
    tt[j] <- tau
    ev[j] <- 0L
  }
  data.frame(time = tt, event = ev, censored = 1L - ev, trt = trt,
             centre = factor(centre))
}

for (nd in c(5L, 0L)) {
  cat(sprintf("\n== n_nodeath %d ==\n", nd))
  for (s in 20260910L + 0:5) {
    d <- sim_nodeath(s, n_nodeath = nd)
    bk <- range(log(d$time[d$event == 1L]))
    f <- suppressWarnings(frm(
      bf(time | cens(censored) ~ trt, gamma1 ~ (1 | centre)),
      family = royston_parmar(knots = numeric(0), bknots = bk),
      data = d, se = TRUE))
    g1 <- unname(frmtmb::fixef_by_dpar(f)$gamma1[["(Intercept)"]])
    sh <- g1 + as.numeric(frmtmb::ranef(f)[["centre"]])
    # the warning at fit end, captured by refitting with the handler on
    w <- character(0)
    withCallingHandlers(
      frm(bf(time | cens(censored) ~ trt, gamma1 ~ (1 | centre)),
          family = royston_parmar(knots = numeric(0), bknots = bk),
          data = d, se = TRUE),
      warning = function(cnd) {
        w <<- c(w, conditionMessage(cnd))
        invokeRestart("muffleWarning")
      })
    r <- rp_floored(f, action = "report")
    ref <- tryCatch({rp_floored(f); "no"}, error = function(e) "yes")
    cv <- tryCatch({
      frm_curve(f, newdata = data.frame(
        time = exp(seq(-2, 1, length.out = 5)), trt = 0,
        centre = factor(1, levels = levels(d$centre))),
        dpar = "gamma1", simultaneous = FALSE)
      "no"
    }, error = function(e) "yes")
    ncs <- if (is.null(r[["n_nonmonotone_censored"]])) NA_integer_ else
      r[["n_nonmonotone_censored"]]
    cat(sprintf(paste0("seed %d  min slope %8.4f  centres<=0 %d  ",
                       "n_nonmonotone %d  n_nonmonotone_censored %s  ",
                       "refuses %s  frm_curve refuses %s  ",
                       "fit-end warnings %d\n"),
                s, min(sh), sum(sh <= 0), r[["n_nonmonotone"]],
                format(ncs), ref, cv, length(w)))
    rr <- attr(r, "rows")[["nonmonotone_censored"]]
    nd_rows <- which(as.integer(d$centre) <= nd)
    neg <- sh[sh <= 0]
    cat("  slopes <= 0:",
        if (length(neg)) sprintf("%.4f", sort(neg)) else "none",
        " rows flagged", if (is.null(rr)) "NA" else length(rr),
        " equal to the no-death rows",
        if (is.null(rr)) "NA" else identical(as.integer(rr), nd_rows),
        " max abs gradient", sprintf("%.2e", frmtmb::diagnose(f, quiet = TRUE)$max_grad),
        " convergence", frmtmb::diagnose(f, quiet = TRUE)$convergence, "\n")
    if (length(w)) cat("  warning:", substr(w[1L], 1L, 120L), "\n")
  }
}
