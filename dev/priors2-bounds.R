# Lane priors2: prior bounds on the fit route (backlog 5005, 5019, 5027)
# and a bound on class Intercept (5646), before and after.
# PRIORS2_LIB unset = lane library; PRIORS2_LIB=none = rellib-r21 only.
# Rscript dev/priors2-bounds.R > dev/priors2-out/bounds-<arm>.txt
lib <- Sys.getenv("PRIORS2_LIB", "C:/Users/adf44/source/r/wt-priors2-lib")
if (identical(lib, "none")) lib <- ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
options(width = 150, digits = 10)
show <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", substr(gsub("\n", " ", conditionMessage(w)), 1, 160),
               "\n")
           invokeRestart("muffleWarning")
         },
         message = function(m) invokeRestart("muffleMessage")),
       error = function(e) {
         cat("ERROR: ", substr(gsub("\n", " ", conditionMessage(e)), 1, 300),
             "\n")
         NULL
       })
  if (!is.null(r)) print(r)
  invisible(r)
}
summ <- function(f, nm = NULL) {
  out <- c(conv = f$opt$convergence, obj = f$opt$objective)
  if (!is.null(nm)) out <- c(out, f$opt$par[nm])
  out
}
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$yneg <- d$y - 0.8 * d$z
d$ys <- 1 + rnorm(n, sd = exp(0.2 + 0.3 * d$x))
zf <- function(f) unname(fixef(f)["z", "Estimate"])

## 5005: a user's lb where the density is 0 or infinite
for (dn in c("lognormal(0, 1)", "beta(2, 2)", "gamma(0.5, 1)",
             "gamma(2, 2)", "exponential(1)")) {
  show(paste("A", dn, "lb = 0 on b z; yneg (data want z < 0)"), {
    f <- frm(yneg ~ z, data = d,
             prior = set_prior(dn, class = "b", coef = "z", lb = 0,
                               ub = if (grepl("beta", dn)) 1 else NA))
    c(summ(f), z = zf(f))
  })
  show(paste("A", dn, "lb = 1e-8 on b z; yneg"), {
    f <- frm(yneg ~ z, data = d,
             prior = set_prior(dn, class = "b", coef = "z", lb = 1e-8,
                               ub = if (grepl("beta", dn)) 1 else NA))
    c(summ(f), z = zf(f))
  })
  show(paste("A", dn, "lb = 0 on b z; y (interior mode)"), {
    f <- frm(y ~ z, data = d,
             prior = set_prior(dn, class = "b", coef = "z", lb = 0,
                               ub = if (grepl("beta", dn)) 1 else NA))
    c(summ(f), z = zf(f))
  })
}

## 5019: uniform on a centered intercept
show("B7 uniform(0, 0.1) on sigma's Intercept (sigma ~ x)", {
  f <- frm(bf(ys ~ 1, sigma ~ x), data = d,
           prior = set_prior("uniform(0, 0.1)", class = "Intercept",
                             dpar = "sigma"))
  summ(f)
})
show("B10 uniform(0, 0.5) on centered mu Intercept", {
  f <- frm(y ~ z, data = d, prior = set_prior("uniform(0, 0.5)",
                                             class = "Intercept"))
  e <- fixef(f)[, "Estimate"]
  c(summ(f), e, centered = unname(e[1] + e[2] * mean(d$z)))
})
show("B10b uniform(0, 0.5) on centered mu Intercept, lb/ub given", {
  f <- frm(y ~ z, data = d, prior = set_prior("uniform(0, 0.5)",
                                             class = "Intercept",
                                             lb = 0, ub = 0.5))
  e <- fixef(f)[, "Estimate"]
  c(summ(f), e, centered = unname(e[1] + e[2] * mean(d$z)))
})

## 5027: a box on a dpar's own class one ulp outside the support
set.seed(1)
dz <- data.frame(z = rnorm(n))
invisible(rnorm(10)); invisible(rnorm(n))
dz$cnt <- rpois(n, 2) * rbinom(n, 1, 0.5)
show("C1 zi ML", summ(frm(cnt ~ z, data = dz,
                          family = zero_inflated_poisson())))
show("C2 zi uniform(0, 0.3), no bound", {
  f <- frm(cnt ~ z, data = dz, family = zero_inflated_poisson(),
           prior = set_prior("uniform(0, 0.3)", class = "zi"))
  c(summ(f), zi = plogis(unname(f$estimates$betad[1])))
})
show("C3 zi uniform(0, 0.3), lb 0 ub 0.3", {
  f <- frm(cnt ~ z, data = dz, family = zero_inflated_poisson(),
           prior = set_prior("uniform(0, 0.3)", class = "zi", lb = 0,
                             ub = 0.3))
  c(summ(f), zi = plogis(unname(f$estimates$betad[1])))
})
w <- c(seq(0.31, 0.99, by = 0.01), seq(1.1, 3, by = 0.1))
below <- w[exp(log(w)) < w]
cat("\nlb with exp(log(lb)) < lb:", head(below, 5), "\n")
show("C4 sigma uniform(lb, 5), lb with exp(log(lb)) < lb", {
  lb <- below[1]
  f <- frm(y ~ z, data = d,
           prior = set_prior(sprintf("uniform(%.17g, 5)", lb),
                             class = "sigma"))
  c(summ(f), sigma = exp(unname(f$estimates$betad[1])), lb = lb)
})
show("C5 the same with lb/ub given", {
  lb <- below[1]
  f <- frm(y ~ z, data = d,
           prior = set_prior(sprintf("uniform(%.17g, 5)", lb),
                             class = "sigma", lb = lb, ub = 5))
  c(summ(f), sigma = exp(unname(f$estimates$betad[1])), lb = lb)
})
set.seed(1)
d3 <- data.frame(z = rnorm(n))
d3$y <- 1 + 0.3 * d3$z + rnorm(n, sd = 0.3)
for (lbs in c(0.35, 0.36)) {
  show(paste("C6 sigma uniform(lb, 5), lb", lbs, ", data's sigma below it"), {
    f <- frm(y ~ z, data = d3,
             prior = set_prior(sprintf("uniform(%.17g, 5)", lbs),
                               class = "sigma", lb = lbs, ub = 5))
    s <- exp(unname(f$estimates$betad[1]))
    c(summ(f), sigma = s, inside = s >= lbs)
  })
}
# how often the map lands outside, per link, over 100 ends
ends <- seq(0.01, 0.99, by = 0.01)
cat("\nlogit: plogis(qlogis(v)) > v for", sum(plogis(qlogis(ends)) > ends),
    "of", length(ends), "; < v for", sum(plogis(qlogis(ends)) < ends), "\n")
ends2 <- seq(0.05, 5, by = 0.05)
cat("log: exp(log(v)) < v for", sum(exp(log(ends2)) < ends2), "of",
    length(ends2), "; > v for", sum(exp(log(ends2)) > ends2), "\n")
cat("\ndone\n")
