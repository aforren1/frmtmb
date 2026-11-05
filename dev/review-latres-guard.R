# Attack 7: guards with the condition absent. A nonlinear dpar beside a linear
# mu, and the multi-category formulas that reach mu without a suffix.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
set.seed(20261102)
d <- expand.grid(t = 1:5, g = factor(1:10))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
d$y <- rnbinom(n, size = 3, mu = exp(1 + 0.3 * d$x))
d$un <- rbeta(n, 2, 2)
d$a2 <- rbeta(n, 2, 2) * (1 - d$un)
d$a3 <- 1 - d$un - d$a2
d$bi <- rbinom(n, 5, .4)
d$bj <- rbinom(n, 5, .4)
d$tr <- 10L
d$gs <- rnorm(n)
shw <- function(lab, e) {
  r <- msg_of(eval.parent(substitute(e)))
  cat("---", lab, ":", if (!is.null(r$err)) paste("ERR:", substr(r$err, 1, 260)) else "ok", "\n")
  if (length(r$warn)) cat("  WARN:", substr(r$warn[1], 1, 200), "\n")
  r$val
}
f1 <- shw("nlf(shape), ar on mu", frm(bf(y ~ x + ar(t, g, cov = TRUE)) + nlf(shape ~ exp(a)) + lf(a ~ 1), family = negbinomial(), data = d))
if (!is.null(f1)) print(summary(f1))
f2 <- shw("zi nlf, ar on mu", frm(bf(y ~ x + ar(t, g, cov = TRUE)) + nlf(zi ~ inv_logit(a)) + lf(a ~ 1), family = zero_inflated_poisson(), data = d))
f3 <- shw("no ar, nlf(shape)", frm(bf(y ~ x) + nlf(shape ~ exp(a)) + lf(a ~ 1), family = negbinomial(), data = d))
for (spec in list(
  list("dir 2 cats", bf(cbind(un, oth) ~ x + ar(t, g, cov = TRUE)), dirichlet()),
  list("dir 3 cats", bf(cbind(un, a2, a3) ~ x + ar(t, g, cov = TRUE)), dirichlet()),
  list("dir2 3 cats", bf(cbind(un, a2, a3) ~ x + ar(t, g, cov = TRUE)), "dirichlet2"),
  list("logistic_normal 3", bf(cbind(un, a2, a3) ~ x + ar(t, g, cov = TRUE)), frmtmb::logistic_normal()),
  list("dirmult", bf(cbind(bi, bj) | trials(tr) ~ x + ar(t, g, cov = TRUE)), dirichlet_multinomial()),
  list("dir explicit mu", bf(cbind(un, oth) ~ 1, muun ~ x + ar(t, g, cov = TRUE)), dirichlet())
)) {
  dd <- d
  dd$oth <- 1 - dd$un
  r <- msg_of(frm(spec[[2]], family = spec[[3]], data = dd, dry_run = "frame"))
  cat("---", spec[[1]], ":", if (!is.null(r$err)) paste("ERR:", substr(r$err, 1, 200)) else "BUILDS", "\n")
  if (is.null(r$err)) {
    fr <- r$val
    cat("   autocor:", length(fr$autocor), " latres:", length(fi("frame_latres")(fr)), " re_blocks:", length(fr$re_blocks), "\n")
    print(names(fr$linpreds))
  }
}
cat("\n== the nl guard's wording under nl = TRUE autocor spelling\n")
shw("nl autocor poisson", frm(bf(y ~ a * exp(b * x), a + b ~ 1, nl = TRUE, autocor = ~ ar(t, g, cov = TRUE)), family = poisson(), data = d))
shw("autocor= arg, linear", frm(bf(y ~ x, autocor = ~ ar(t, g, cov = TRUE)), family = poisson(), data = d))
