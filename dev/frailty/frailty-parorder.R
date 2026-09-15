# The parameter-ordering guard, run against the case it exists to catch.
#
# `frailty-sweep.R` and `frailty-ident.R` score the Laplace objective at
# a named parameter point through
# `setNames(c(gam[1], beta, gam[-1], log(sd)), names(dry$obj$par))`.
# `setNames()` relabels and does not reorder, so that line ASSERTS an
# ordering, and the objective's own names are `beta, beta, betad, betad,
# theta`: duplicated, so they cannot protect it. Every offset in
# dev/frailty-findings.md is scored through that line.
#
# Raised by the punch round; the review found the assumption right and
# unasserted (dev/reviews/20260910-frailty.md, section 2e).
#
# SEED 20260910.
source("frailty-common.R")
suppressMessages({
  library(frmtmb); library(frmtmb.spline)
})

d <- frailty_sim(20260910L)
kn <- rp_knots(log(d$time[d$event == 1L]), 2L)
fam <- royston_parmar(knots = kn$ik, bknots = kn$bk)
form <- bf(time | cens(censored) ~ trt + (1 | centre))
ff <- frm(form, family = fam, data = d, se = TRUE)
dry <- frm(form, family = fam, data = d, dry_run = "objective")

gam <- c(unname(fixef(ff)$mu[["(Intercept)"]]),
         unname(fixef(ff)$gamma1), unname(fixef(ff)$gamma2))
beta <- unname(fixef(ff)$mu[["trt"]])
sd_b <- sqrt(frmtmb::VarCorr(ff)[[1L]][1L, 1L])
ll <- as.numeric(logLik(ff))

cat("obj$par names:", paste(names(dry$obj$par), collapse = ", "), "\n")
cat("duplicated:", any(duplicated(names(dry$obj$par))), "\n\n")

score <- function(p) -dry$obj$fn(stats::setNames(p, names(dry$obj$par)))

# as the scripts build it
right <- c(gam[1L], beta, gam[-1L], log(sd_b))
# every other ordering of the same five numbers that a plausible slip
# would produce: the covariate ahead of the intercept, the gammas
# swapped, and the two blocks exchanged
wrong <- list(
  `beta ahead of the intercept` = c(beta, gam[1L], gam[-1L], log(sd_b)),
  `the two gammas swapped` = c(gam[1L], beta, gam[3L], gam[2L],
                               log(sd_b)),
  `mu block and gamma block exchanged` = c(gam[-1L], gam[1L], beta,
                                           log(sd_b)))

cat(sprintf("%-38s %16s %12s\n", "ordering", "-obj$fn", "rel to logLik"))
cat(sprintf("%-38s %16.8f %12.2e   <- as shipped\n", "as the scripts build it",
            score(right), abs(score(right) - ll) / abs(ll)))
for (nm in names(wrong)) {
  v <- score(wrong[[nm]])
  cat(sprintf("%-38s %16.8f %12.2e\n", nm, v, abs(v - ll) / abs(ll)))
}
cat("\nlogLik(ff) =", format(ll, digits = 12), "\n")
cat("the guard's threshold is 1e-8 relative; every wrong ordering is\n",
    "orders of magnitude outside it, so the guard fails CLOSED\n",
    sep = "")
