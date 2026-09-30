# Reviewer, claims 2 to 4 against brms 2.23.0's own compiled programs
# (brms_lp_check of tests/testthat/helper-brms.R: brms's log_prob at
# frmtmb's estimates; joint = TRUE is check C, the inner block).
# Seeds inside each case. Cache: the session scratchpad.
# Log: dev/aterms2-rev-log-03-brms.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-aterms2"
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = "C:/Users/adf44/AppData/Local/Temp/1/claude/c--Users-adf44-source-r-frmtmb/66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad/aterms2-rev-stan-cache",
           FRMTMB_BRMS_FIT_TESTS = "true", NOT_CRAN = "true")

suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
source(file.path(wt, "tests/testthat/helper-brms.R"))
only <- commandArgs(TRUE)
run <- function(label, expr) {
  if (length(only) && !label %in% only) return(invisible(NULL))
  cat("\n==========", label, "\n")
  r <- tryCatch(expr, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    NULL
  }, expectation_failure = function(e) {
    cat("EXPECTATION FAILED:", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(r)) {
    cat(sprintf("measured_const %.6g max_grad %.3g ours %.10g\n",
                r$measured_const, r$max_grad, r$ours))
  }
  invisible(r)
}
q <- function(expr) suppressWarnings(suppressMessages(expr))

## ---- claim 2: |ID| shared across a subsetted and an unsubsetted response
set.seed(301)
n <- 160
G <- 8
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(rep(seq_len(G), length.out = n)))
L <- t(chol(matrix(c(1, 0.6, 0.6, 1), 2)))
u <- t(L %*% matrix(rnorm(2 * G), 2)) * 0.7
d$y1 <- 1 + d$x + u[d$g, 1] + rnorm(n)
d$y2 <- rpois(n, exp(0.3 - 0.4 * d$z + u[d$g, 2]))
d$s1 <- d$x > -0.4
stopifnot(all(table(d$g[d$s1]) > 0))
d$s2 <- rep(c(TRUE, TRUE, FALSE), length.out = n)
stopifnot(all(table(d$g[d$s2]) > 0))
fid <- q(frm(bf(y1 | subset(s1) ~ x + (1 | p | g)) + gaussian() +
               bf(y2 ~ z + (1 | p | g)) + poisson(), data = d))
bid <- brms::bf(y1 | subset(s1) ~ x + (1 | p | g), family = gaussian()) +
  brms::bf(y2 ~ z + (1 | p | g), family = poisson()) +
  brms::set_rescor(FALSE)
run("id_subset_one", brms_lp_check(bid, NULL, d, fid, joint = TRUE))
fid2 <- q(frm(bf(y1 | subset(s1) ~ x + (1 + x | p | g)) + gaussian() +
                bf(y2 | subset(s2) ~ z + (1 | p | g)) + poisson(), data = d))
bid2 <- brms::bf(y1 | subset(s1) ~ x + (1 + x | p | g),
                 family = gaussian()) +
  brms::bf(y2 | subset(s2) ~ z + (1 | p | g), family = poisson()) +
  brms::set_rescor(FALSE)
run("id_subset_both", brms_lp_check(bid2, NULL, d, fid2, joint = TRUE))

## ---- claim 2: interval censoring (cens_y2 rows) with subset and NA drop
set.seed(302)
n <- 120
dc <- data.frame(x = rnorm(n), z = rnorm(n),
                 s1 = rep(c(TRUE, FALSE, TRUE), length.out = n),
                 s2 = rep(c(TRUE, TRUE, FALSE, TRUE), length.out = n))
dc$y1 <- 1 + dc$x + rnorm(n)
dc$cc <- rep(c(0, 2, 1, 0, -1), length.out = n)
dc$yup <- dc$y1 + runif(n, 0.3, 1)
dc$y2 <- 0.5 * dc$z + rnorm(n)
# dropped for every response: x is NA on rows 4 and 10, both inside s1
dc$x[c(4, 10)] <- NA
# harmless: y2 and z NA only where s2 is FALSE
dc$y2[which(!dc$s2)[1:5]] <- NA
dc$z[which(!dc$s2)[6:8]] <- NA
fc <- q(frm(bf(y1 | subset(s1) + cens(cc, yup) ~ x) +
              bf(y2 | subset(s2) ~ z), data = dc, family = gaussian()))
cat("cens model n_obs", fc$frame$n_obs, " rows y1", length(fc$frame$y$y1),
    " rows y2", length(fc$frame$y$y2), "\n")
bc <- brms::bf(y1 | subset(s1) + cens(cc, yup) ~ x) +
  brms::bf(y2 | subset(s2) ~ z) + brms::set_rescor(FALSE)
run("subset_cens_interval_na", brms_lp_check(bc, gaussian(), dc, fc))

## ---- claim 2: subset with weights, trials and trunc in one model
set.seed(303)
n <- 150
dw <- data.frame(x = rnorm(n), z = rnorm(n), wt = runif(n, 0.5, 2),
                 nt = sample(5:12, n, TRUE),
                 s1 = rep(c(TRUE, FALSE), length.out = n),
                 s2 = rep(c(FALSE, TRUE, TRUE), length.out = n))
dw$y1 <- pmax(1 + dw$x + rnorm(n), -0.5)
dw$yb <- rbinom(n, dw$nt, plogis(0.2 + 0.6 * dw$z))
fw <- q(frm(bf(y1 | subset(s1) + weights(wt) + trunc(lb = -0.5) ~ x) +
              gaussian() + bf(yb | subset(s2) + trials(nt) ~ z) +
              binomial(), data = dw))
bw <- brms::bf(y1 | subset(s1) + weights(wt) + trunc(lb = -0.5) ~ x,
               family = gaussian()) +
  brms::bf(yb | subset(s2) + trials(nt) ~ z, family = binomial()) +
  brms::set_rescor(FALSE)
run("subset_weights_trials_trunc", brms_lp_check(bw, NULL, dw, fw))

## ---- claim 2: smooth, mo() and a distributional sigma under subset
set.seed(304)
n <- 160
ds <- data.frame(x = runif(n, -2, 2), z = rnorm(n),
                 om = factor(sample(1:4, n, TRUE), ordered = TRUE),
                 s1 = rep(c(TRUE, TRUE, FALSE), length.out = n))
ds$y1 <- sin(ds$x) + as.integer(ds$om) * 0.3 + rnorm(n, 0, exp(0.2 * ds$z))
ds$y2 <- ds$z + rnorm(n)
fsm <- q(frm(bf(y1 | subset(s1) ~ mo(om) + x, sigma ~ z) + bf(y2 ~ z),
             data = ds, family = gaussian()))
bsm <- brms::bf(y1 | subset(s1) ~ mo(om) + x, sigma ~ z) +
  brms::bf(y2 ~ z) + brms::set_rescor(FALSE)
run("subset_mo_sigma", brms_lp_check(bsm, gaussian(), ds, fsm))

## ---- claim 3: mi(x, idx = ) with the rows of x in a different order,
## character ids, and x subsetted
set.seed(305)
n <- 200
ids <- paste0("id", sprintf("%03d", sample(1000, n)))
dm <- data.frame(xid = ids, s = rep(c(TRUE, FALSE, TRUE, TRUE),
                                    length.out = n), w = rnorm(n))
dm$x <- rnorm(n)
# y's rows point at x's rows in a shuffled order, many-to-one
xrows <- which(dm$s)
dm$ref <- dm$xid[sample(xrows, n, TRUE)]
dm$y <- 1 + 0.6 * dm$x[match(dm$ref, dm$xid)] + 0.3 * dm$w +
  rnorm(n, sd = 0.5)
dm$x[xrows[c(2, 11, 30)]] <- NA
fm <- q(frm(bf(y ~ mi(x, idx = ref) + w) +
              bf(x | mi() + index(xid) + subset(s) ~ 1), data = dm,
            family = gaussian()))
bm <- brms::bf(y ~ mi(x, idx = ref) + w) +
  brms::bf(x | mi() + index(xid) + subset(s) ~ 1) + brms::set_rescor(FALSE)
sdb <- brms::standata(bm, data = dm)
mt <- fm$frame$linpreds[[which(vapply(fm$frame$linpreds, function(l)
  length(l$mi) > 0, TRUE))]]$mi[[1]]
cat("idxl identical to brms's idxl_y_x_1:",
    identical(as.integer(mt$idxl), as.integer(sdb$idxl_y_x_1)),
    " length", length(mt$idxl), "\n")
cat("brms Jmi_x:", sdb$Jmi_x, "\n")
run("mi_idx_shuffled", brms_lp_check(bm, gaussian(), dm, fm, joint = TRUE))

## ---- claim 4: rate() under the links and terms the worker did not run
set.seed(306)
n <- 250
dr <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4),
                 wt = runif(n, 0.5, 2))
dr$y <- rpois(n, exp(0.4 + 0.3 * dr$x) * dr$time)
dr$yn <- rnbinom(n, mu = (2 + 0.5 * dr$x^2) * dr$time, size = 3 * dr$time)
dr$cc <- rep(c(0, 0, 1, -1), length.out = n)
dr$yt <- pmax(dr$y, 1L)
fr1 <- q(frm(yn | rate(time) ~ 1 + I(x^2), data = dr,
             family = negbinomial("identity")))
run("rate_negbin_identity",
    brms_lp_check(brms::bf(yn | rate(time) ~ 1 + I(x^2)),
                  brms::negbinomial("identity"), dr, fr1))
fr2 <- q(frm(yn | rate(time) ~ 1 + I(x^2), data = dr,
             family = geometric("identity")))
run("rate_geometric_identity",
    brms_lp_check(brms::bf(yn | rate(time) ~ 1 + I(x^2)),
                  brms::geometric("identity"), dr, fr2))
fr3 <- q(frm(yn | rate(time) ~ x, data = dr, family = negbinomial("sqrt")))
run("rate_negbin_sqrt",
    brms_lp_check(brms::bf(yn | rate(time) ~ x),
                  brms::negbinomial("sqrt"), dr, fr3))
fr4 <- q(frm(y | rate(time) + cens(cc) ~ x, data = dr, family = poisson()))
run("rate_poisson_cens",
    brms_lp_check(brms::bf(y | rate(time) + cens(cc) ~ x), poisson(), dr,
                  fr4))
fr5 <- q(frm(yt | rate(time) + trunc(lb = 1) ~ x, data = dr,
             family = poisson()))
run("rate_poisson_trunc",
    brms_lp_check(brms::bf(yt | rate(time) + trunc(lb = 1) ~ x), poisson(),
                  dr, fr5))
fr6 <- q(frm(y | rate(time) + weights(wt) ~ x, data = dr,
             family = poisson()))
run("rate_poisson_weights",
    brms_lp_check(brms::bf(y | rate(time) + weights(wt) ~ x), poisson(),
                  dr, fr6))
fr7 <- q(frm(yn | rate(time) ~ x, data = dr, family = geometric()))
run("rate_geometric_log_again",
    brms_lp_check(brms::bf(yn | rate(time) ~ x), brms::geometric(), dr,
                  fr7))
