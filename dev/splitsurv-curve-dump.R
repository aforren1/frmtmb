# The output of frm_curve(), frm_curve_deriv() and frm_curve_feature()
# on the RP and non-RP cases of the moved tests, saved so that the base
# build (frmtmb.spline, rellib-r11) and the lane's build (core plus
# frmtmb.survival) can be compared value for value.
#
#   Rscript dev/splitsurv-curve-dump.R <before|after> <out.rds>
#
# "before" reads rellib-r11 alone; "after" puts the lane's library
# first. Each case runs in this one process, so a difference between
# the two arms is a difference between builds and not between cores
# (dev/rtmb-pitfalls.md, item 21, is answered by repeating the pair).
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
out <- a[2]
LIB <- "C:/Users/adf44/source/r/wt-splitsurv-lib"
REL <- "C:/Users/adf44/source/r/rellib-r11"
USR <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(if (identical(arm, "after")) LIB, REL, USR))
if (identical(arm, "after")) {
  suppressMessages(library(frmtmb.survival))
} else {
  suppressMessages(library(frmtmb.spline))
}
cat("frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("frm_curve from", environmentName(environment(frm_curve)), "\n")

res <- list()
keep <- function(nm, expr) {
  v <- tryCatch(expr, error = function(e) {
    structure(conditionMessage(e), class = "dump_error")
  })
  res[[nm]] <<- v
  invisible(v)
}

# ---- non-RP: test-curve.R's sp_curve_fit()
set.seed(7)
d <- data.frame(x = sort(stats::runif(250)))
d$y <- 2 * sin(pi * d$x) + 0.6 * d$x + stats::rnorm(250, 0, 0.4)
fit <- frm(bf(y ~ s(x, k = 10)), family = gaussian(), data = d)
g <- data.frame(x = seq(0, 1, length.out = 40))
keep("curve", frm_curve(fit, newdata = g, nsim = 2000, seed = 1))
keep("curve_transform", frm_curve(fit, newdata = g, simultaneous = FALSE,
                                  transform = TRUE))
keep("deriv1", frm_curve_deriv(fit, var = "x", order = 1, newdata = g,
                               nsim = 2000, seed = 2))
keep("deriv2", frm_curve_deriv(fit, var = "x", order = 2, newdata = g,
                               simultaneous = FALSE))
gf <- data.frame(x = seq(0.05, 0.95, length.out = 41))
keep("feat_max", frm_curve_feature(fit, var = "x", type = "maximum",
                                   newdata = gf))
keep("feat_cross", frm_curve_feature(fit, var = "x", type = "crossing",
                                     at = 1, newdata = gf[1:15, , drop = FALSE]))

# ---- non-RP: a difference curve (test-difference.R's design)
set.seed(3)
n <- 400
dd <- data.frame(x = stats::runif(n), fac = factor(rep(c("a", "b"), n / 2)))
dd$y <- sin(2 * pi * dd$x) + (dd$fac == "b") * dd$x + stats::rnorm(n, 0, 0.3)
fd <- frm(bf(y ~ fac + s(x, by = fac, k = 8)), family = gaussian(),
          data = dd)
ga <- data.frame(x = seq(0, 1, length.out = 30),
                 fac = factor("a", levels = levels(dd$fac)))
gb <- data.frame(x = seq(0, 1, length.out = 30),
                 fac = factor("b", levels = levels(dd$fac)))
keep("contrast", frm_curve(fd, newdata = gb, contrast = ga, nsim = 2000,
                           seed = 3))
keep("contrast_deriv", frm_curve_deriv(fd, var = "x", newdata = gb,
                                       contrast = ga, simultaneous = FALSE))

# ---- RP: test-rp-floored.R's sp_far_censored()
set.seed(20260905)
dr <- data.frame(grp = factor(rep(c("A", "B"), each = 300)))
dr$t <- stats::rweibull(600, shape = 1.6,
                        scale = ifelse(dr$grp == "A", 0.30, 3.0))
dr$censored <- 0L
i <- which(dr$grp == "A")[1L]
dr$t[i] <- 50
dr$censored[i] <- 1L
fr <- suppressWarnings(frm(bf(t | cens(censored) ~ grp),
                           family = royston_parmar(df = 3), data = dr))
gr <- data.frame(t = seq(0.05, 5, length.out = 10),
                 grp = factor("A", levels = levels(dr$grp)))
keep("rp_curve_mu", frm_curve(fr, newdata = gr, dpar = "mu", nsim = 2000,
                              seed = 4))
keep("rp_deriv_mu", frm_curve_deriv(fr, var = "t", newdata = gr,
                                    dpar = "mu", simultaneous = FALSE))
keep("rp_floored_report", rp_floored(fr, action = "report"))
# the gate: a non-monotone fit refuses every curve function
bad <- fr
bad$estimates$betad[] <- bad$estimates$betad - 3
keep("rp_bad_report", rp_floored(bad, action = "report"))
keep("rp_bad_curve", frm_curve(bad, newdata = gr, dpar = "mu",
                               simultaneous = FALSE))
keep("rp_bad_deriv", frm_curve_deriv(bad, var = "t", newdata = gr,
                                     dpar = "mu", simultaneous = FALSE))
keep("rp_bad_feature", frm_curve_feature(bad, var = "t", newdata = gr,
                                         dpar = "mu", type = "maximum"))

# ---- RP: flexsurv's bc, the gamma1 curve and a feature on it
e <- new.env()
utils::data("bc", package = "flexsurv", envir = e)
bc <- e$bc
bc$censored <- 1 - bc$censrec
fb <- frm(bf(recyrs | cens(censored) ~ group),
          family = royston_parmar(df = 2), data = bc)
gbc <- data.frame(recyrs = seq(0.5, 6, length.out = 12),
                  group = factor("Good", levels = levels(bc$group)))
keep("bc_curve_gamma1", frm_curve(fb, newdata = gbc, dpar = "gamma1",
                                  simultaneous = FALSE))
keep("bc_curve_mu", frm_curve(fb, newdata = gbc, dpar = "mu", nsim = 2000,
                              seed = 5))
keep("bc_deriv_mu", frm_curve_deriv(fb, var = "recyrs", newdata = gbc,
                                    dpar = "mu", nsim = 2000, seed = 6))
keep("bc_deriv_mu_pw", frm_curve_deriv(fb, var = "recyrs", newdata = gbc,
                                       dpar = "mu", simultaneous = FALSE))

saveRDS(res, out)
for (nm in names(res)) {
  v <- res[[nm]]
  cat(sprintf("%-18s %s\n", nm, if (inherits(v, "dump_error")) {
    paste("ERROR:", substr(unclass(v), 1, 90))
  } else class(v)[1]))
}
