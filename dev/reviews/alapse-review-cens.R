# Review of lane alapse, check 3: ald_cens_sep() guard-absent cases
# beyond the lane's 72 of 80. Data of dev/wcheck2-allcens.R (seed
# 20261021, 80 rows, right censored at their minimum).
#   Rscript dev/reviews/alapse-review-cens.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/reviews/alapse-review-common.R")
set.seed(20261021)
n <- 80
d0 <- data.frame(x = stats::rnorm(n), g4 = factor(rep(1:20, each = 4)))
d0$ya <- 1 + 0.5 * d0$x + stats::rnorm(n)
d0$ca <- "right"
d0$yac <- rep(min(d0$ya), n)
d0$w <- 1

run <- function(id, fo, d, fam = asym_laplace(), prior = NULL) {
  cat("==", id, "\n")
  r <- cap(frm(fo, family = fam, data = d, prior = prior))
  say("frm", r)
  if (inherits(r$v, "error")) return(invisible())
  f <- r$v
  sep <- if (identical(lib_arg, "lane")) ns$ald_cens_sep(f)
  cat(sprintf("  code %d | logLik %.6g | ald_cens_sep %s | sep warn %s\n",
              f$opt$convergence, as.numeric(logLik(f)),
              if (is.null(sep)) "NULL" else "fires",
              "frmtmb_separation" %in% r$wc))
  print(round(f$opt$par, 3))
  lost <- ns$sdr_of(f)$se_lost
  cat("  lost:", paste(names(lost), lost, sep = "=", collapse = " "), "\n")
}

run("all-right", yac | cens(ca) ~ x, d0)
d <- d0; d$ca <- "left"
d$yac <- rep(max(d0$ya), n)
run("all-left", yac | cens(ca) ~ x, d)
# one exact row, below and above the censoring point
d <- d0; d$ca[1] <- "none"; d$yac[1] <- min(d0$ya) - 1
run("one-exact-below", yac | cens(ca) ~ x, d)
d <- d0; d$ca[1] <- "none"; d$yac[1] <- min(d0$ya) + 1
run("one-exact-above", yac | cens(ca) ~ x, d)
# one exact row with weight 0: the proof still holds
d <- d0; d$ca[1] <- "none"; d$yac[1] <- min(d0$ya) - 1; d$w[1] <- 0
run("one-exact-w0", yac | cens(ca) + weights(w) ~ x, d)
# censored both sides
d <- d0; d$ca[41:80] <- "left"; d$yac[41:80] <- max(d0$ya)
run("both-sides-apart", yac | cens(ca) ~ x, d)
d <- d0; d$ca[41:80] <- "left"; d$yac[1:40] <- d0$ya[1:40] + 0.5
d$yac[41:80] <- d0$ya[41:80] - 0.5
run("both-sides-overlap", yac | cens(ca) ~ x, d)
# no intercept; a slope only
run("no-intercept", yac | cens(ca) ~ 0 + x, d0)
# priors: on the intercept, and on the slope only
run("prior-intercept", yac | cens(ca) ~ x, d0,
    prior = set_prior("normal(0, 1)", class = "Intercept"))
run("prior-slope", yac | cens(ca) ~ x, d0,
    prior = set_prior("normal(0, 1)", class = "b"))
# fixed quantile; a group term in sigma
run("fixed-quantile", bf(yac | cens(ca) ~ x, quantile = 0.2), d0)
run("sigma-group", bf(yac | cens(ca) ~ x, sigma ~ (1 | g4)), d0)
# interval censoring everywhere
d <- d0; d$ca <- "interval"; d$y2 <- d$yac + 1
run("interval", yac | cens(ca, y2) ~ x, d)
# gaussian with every row censored (the gaussian path is untouched)
run("gaussian-all-right", yac | cens(ca) ~ x, d0, fam = gaussian())
cat("DONE\n")
