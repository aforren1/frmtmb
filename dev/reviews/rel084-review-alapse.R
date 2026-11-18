# rel084 release review, alapse fix round 2: (1) a finite bound on sigma's
# or the quantile's coefficient disables the search and the fit stays in
# its box; a bound on a location coefficient leaves the search; (2) the
# search restarts through any positive / unit link; (3) the local warning's
# advice follows the reason. New data (skewed errors, seed 8411), not the
# lane's design B.
#   REVLIB=r23 Rscript dev/reviews/rel084-review-alapse.R
source("dev/reviews/rel084-review-common.R")
options(frmtmb.notices = FALSE, width = 160)
ns <- asNamespace("frmtmb")
set.seed(8411)
n <- 120
d <- data.frame(x = rnorm(n))
d$y <- 0.5 + 0.8 * d$x + rexp(n) - 1 + 0.3 * rnorm(n)
inside_box <- function(f) {
  bx <- f$cache$bounds
  np <- length(f$opt$par)
  all(f$opt$par >= rep_len(bx$lower %||% -Inf, np) &
        f$opt$par <= rep_len(bx$upper %||% Inf, np))
}
row <- function(lab, r) {
  f <- r$v
  if (is.null(f)) { cat(sprintf("%-44s ERROR %s\n", lab, cls_line(r))); return(invisible()) }
  dp <- ns$eval_dpars(f)[[1]]
  cat(sprintf("%-44s code %d global %-5s inside %-5s p %.4f sigma %.4f ll %.8f | %s\n",
              lab, f$opt$convergence, isTRUE(f$opt$ald_cert$global),
              inside_box(f), dp$quantile[1], dp$sigma[1],
              as.numeric(logLik(f)), cls_line(r)))
  invisible(f)
}
al <- asym_laplace()
f0 <- row("estimated quantile, no bound",
          conds_of(frm(bf(y ~ x, quantile ~ 1), family = al, data = d)))
q0 <- qlogis(ns$eval_dpars(f0)[[1]]$quantile[1])
s0 <- log(ns$eval_dpars(f0)[[1]]$sigma[1])
cat(sprintf("unbounded optimum logit q %.4f log sigma %.4f\n", q0, s0))
# binding box on the quantile: excludes the unbounded optimum
row("quantile Int box binding [q0-3, q0-0.5]",
    conds_of(frm(bf(y ~ x, quantile ~ 1), family = al, data = d,
                 prior = set_prior("", class = "Intercept", dpar = "quantile",
                                   lb = q0 - 3, ub = q0 - 0.5))))
row("sigma Int box binding [s0+0.3, 3]",
    conds_of(frm(bf(y ~ x, quantile ~ 1, sigma ~ 1), family = al, data = d,
                 prior = set_prior("", class = "Intercept", dpar = "sigma",
                                   lb = s0 + 0.3, ub = 3))))
# a box that does not bind: the proof is off, so the fit should be local
row("quantile Int box wide [-10, 10]",
    conds_of(frm(bf(y ~ x, quantile ~ 1), family = al, data = d,
                 prior = set_prior("", class = "Intercept", dpar = "quantile",
                                   lb = -10, ub = 10))))
# a bound on the slope of mu only
bx <- fixef(f0)
row("mu slope bound [0, 0.3] (binding)",
    conds_of(frm(bf(y ~ x, quantile ~ 1), family = al, data = d,
                 prior = set_prior("", class = "b", coef = "x",
                                   lb = 0, ub = 0.3))))
# links
ll0 <- as.numeric(logLik(f0))
for (lk in c("softplus", "squareplus", "identity")) {
  r <- conds_of(frm(bf(y ~ x, quantile ~ 1),
                    family = asym_laplace(link_sigma = lk), data = d))
  f <- row(paste("sigma link", lk), r)
  if (!is.null(f)) cat(sprintf("   rel diff logLik vs log link %.3e\n",
                               abs(as.numeric(logLik(f)) - ll0) / abs(ll0)))
}
r <- conds_of(frm(bf(y ~ x, quantile ~ 1),
                  family = asym_laplace(link_quantile = "identity"), data = d))
f <- row("quantile link identity", r)
if (!is.null(f)) cat(sprintf("   rel diff logLik vs logit %.3e\n",
                             abs(as.numeric(logLik(f)) - ll0) / abs(ll0)))
# advice
d$y2 <- d$y + rnorm(n)
r <- conds_of(frm(bf(y ~ x) + bf(y2 ~ x) + set_rescor(FALSE),
                  family = asym_laplace(), data = d))
row("two responses", r)
cat("  advice:", substr(gsub("\n", " ", msg_of(r, "frmtmb_ald_local")), 1, 400), "\n")
r <- conds_of(frm(bf(y ~ x, sigma ~ x), family = al, data = d))
row("scale model sigma ~ x", r)
cat("  advice:", substr(gsub("\n", " ", msg_of(r, "frmtmb_ald_local")), 1, 400), "\n")
r <- conds_of(frm(bf(y ~ x, quantile ~ 1), family = al, data = d,
                  prior = set_prior("", class = "Intercept", dpar = "quantile",
                                    lb = q0 - 3, ub = q0 - 0.5)))
cat("  bound advice:", substr(gsub("\n", " ", msg_of(r, "frmtmb_ald_local")), 1, 400), "\n")
cat("DONE\n")
