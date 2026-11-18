# Review of lane alapse, check 5: the sdr_sep_lost() relabel. Separated
# fits, the reasons summary() reads, base against lane.
#   Rscript dev/reviews/alapse-review-relabel.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/reviews/alapse-review-common.R")
cases <- list()
set.seed(4)
dp <- data.frame(x = rnorm(60), z = rnorm(60))
dp$y <- as.integer(dp$x > 0)
cases$bern_complete <- list(fo = y ~ x + z, fam = bernoulli(), d = dp)
set.seed(514)
d5 <- data.frame(x = stats::rnorm(240) * 1e-4, z = stats::rnorm(240))
d5$yb <- as.integer(d5$z > 0)
cases$bern_514 <- list(fo = yb ~ z + x, fam = bernoulli(), d = d5)
cases$bern_514_short <- list(fo = yb ~ z + x, fam = bernoulli(), d = d5,
                             ctl = frmtmb_control(optCtrl = list(eval.max = 40,
                                                                 iter.max = 30),
                                                  restarts = 0))
set.seed(7)
df <- data.frame(f = factor(rep(c("a", "b", "c"), each = 30)),
                 x = rnorm(90))
df$y <- rbinom(90, 1, 0.4)
df$y[df$f == "c"] <- 1
cases$bern_factor <- list(fo = y ~ f + x, fam = bernoulli(), d = df)
cases$probit_factor <- list(fo = y ~ f + x, fam = bernoulli("probit"), d = df)
set.seed(8)
dpz <- data.frame(f = factor(rep(c("a", "b", "c"), each = 30)), x = rnorm(90))
dpz$y <- rpois(90, 2)
dpz$y[dpz$f == "b"] <- 0
cases$pois_zero <- list(fo = y ~ f + x, fam = poisson(), d = dpz)
set.seed(9)
dbn <- data.frame(f = factor(rep(c("a", "b"), each = 20)), x = rnorm(40),
                  n = 5)
dbn$y <- rbinom(40, 5, 0.5)
dbn$y[dbn$f == "b"] <- 5
cases$binom_factor <- list(fo = y | trials(n) ~ f + x, fam = binomial(),
                           d = dbn)
for (nm in names(cases)) {
  cs <- cases[[nm]]
  cat("==", nm, "\n")
  r <- cap(frm(cs$fo, family = cs$fam, data = cs$d,
               control = cs$ctl %||% frmtmb_control()))
  say("frm", r)
  if (inherits(r$v, "error")) next
  f <- r$v
  lost <- ns$sdr_of(f)$se_lost
  cat("  code", f$opt$convergence, "| lost:",
      paste(names(lost), lost, sep = "=", collapse = " "), "\n")
  s <- cap(summary(f))
  cat(paste0("  | ", s$v$se_lost %||% "(none)"), sep = "\n")
  ft <- cap(fitted(f))
  say("fitted", ft)
}
cat("DONE\n")
