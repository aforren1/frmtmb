# rel084 release review, cdfname fix round 1: every reader of the slot
# refuses an old-shaped family (fits other than the lane's tests: censored
# lognormal, truncated negbinomial, uncensored gamma OSA), and `$` returns
# NULL for an absent slot frmtmb_family() defines while it refuses a prefix.
#   REVLIB=r23 Rscript dev/reviews/rel084-review-cdfname.R
source("dev/reviews/rel084-review-common.R")
old_shape <- function(x) {
  if (inherits(x, "frmtmb_family")) {
    names(x)[names(x) == "cdf"] <- "lcdf"
    return(x)
  }
  if (is.list(x) && !is.data.frame(x)) {
    for (i in seq_along(x)) {
      if (is.list(x[[i]])) x[i] <- list(old_shape(x[[i]]))
    }
  }
  x
}
show <- function(lab, r) {
  cat(sprintf("%-46s %s | %s\n", lab, cls_line(r),
              substr(gsub("\n", " ", if (length(r$said))
                r$said[[length(r$said)]]$msg else ""), 1, 90)))
}
set.seed(841)
n <- 150
d <- data.frame(x = rnorm(n))
d$y <- rlnorm(n, 0.2 + 0.5 * d$x, 0.7)
d$cc <- ifelse(d$y > 3, "right", "none")
d$y <- pmin(d$y, 3)
fit_ln <- frm(y | cens(cc) ~ x, family = lognormal(), data = d)
old_ln <- old_shape(fit_ln)
nd <- data.frame(x = c(-1, 0, 1))
show("lognormal cens: predict(newdata) new", conds_of(predict(fit_ln, newdata = nd)))
show("lognormal cens: predict(newdata) old", conds_of(predict(old_ln, newdata = nd)))
show("lognormal cens: update() old", conds_of(update(old_ln, data = d)))
show("lognormal cens: logLik old (no CDF read)", conds_of(logLik(old_ln)))
cat("logLik equal old/new:",
    identical(as.numeric(logLik(old_ln)), as.numeric(logLik(fit_ln))), "\n")
ro <- conds_of(residuals(old_ln, type = "osa"))
rn <- conds_of(residuals(fit_ln, type = "osa"))
show("lognormal cens: osa old", ro)
if (!is.null(ro$v) && !is.null(rn$v)) {
  cat("osa censored old == new:", isTRUE(all.equal(ro$v, rn$v)), "\n")
}

set.seed(842)
dg <- data.frame(x = rnorm(n))
dg$y <- rgamma(n, shape = 3, rate = 3 / exp(0.1 + 0.3 * dg$x))
fit_g <- frm(y ~ x, family = Gamma(link = "log"), data = dg)
old_g <- old_shape(fit_g)
rg <- conds_of(residuals(fit_g, type = "osa"))
cat("gamma osa new finite:", all(is.finite(rg$v[, "Estimate"])), "\n")
show("gamma uncens: osa old", conds_of(residuals(old_g, type = "osa")))
show("gamma uncens: simulate old (no trunc)", conds_of(simulate(old_g, nsim = 2, seed = 1)))

set.seed(843)
dn <- data.frame(x = rnorm(n))
dn$y <- rnbinom(n, size = 3, mu = exp(1 + 0.3 * dn$x))
dn <- dn[dn$y <= 9, ]
fit_nb <- frm(y | trunc(ub = 9) ~ x, family = negbinomial(), data = dn)
old_nb <- old_shape(fit_nb)
show("negbinomial trunc: simulate new", conds_of(simulate(fit_nb, nsim = 2, seed = 1)))
show("negbinomial trunc: simulate old", conds_of(simulate(old_nb, nsim = 2, seed = 1)))

# `$` on the slots
f <- frmtmb_family("my_pois", dpars = "mu", links = list(mu = "log"),
  type = "discrete",
  lpdf = function(y, dpars, aterms) RTMB::dpois(y, dpars[["mu"]], log = TRUE),
  init_dpars = list(mu = function(y, aterms) mean(y)))
cat("custom family has cdf element:", "cdf" %in% names(f), "\n")
f[["cdf"]] <- NULL
r <- conds_of(f$cdf)
cat("absent $cdf:", if (is.null(r$v) && !length(r$said)) "NULL" else cls_line(r), "\n")
r <- conds_of(f$lccdf)
cat("absent/NULL $lccdf:", if (is.null(r$v) && !length(r$said)) "NULL" else cls_line(r), "\n")
show("$cdf_l (prefix)", conds_of(f$cdf_l))
show("$lcdf (old)", conds_of(as_frmtmb_family(poisson())$lcdf))
show("$cd (prefix of cdf and cdf_log)", conds_of(f$cd))
r <- conds_of(as_frmtmb_family(poisson())$zzz)
cat("$zzz:", if (is.null(r$v) && !length(r$said)) "NULL" else cls_line(r), "\n")
show("frmtmb_family(lcdf =)", conds_of(frmtmb_family("x", dpars = "mu",
  links = list(mu = "log"), lpdf = function(y, dpars, aterms) 0,
  lcdf = function(q, dpars, aterms) 0)))
cat("DONE\n")
