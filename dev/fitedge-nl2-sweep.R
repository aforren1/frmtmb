# Item 2 sweep: models with two nonlinear bodies, where the ridge pass of
# nl_fit_off_saddle() can run (a model with one body never reaches it: every
# body is already moved). Data as dev/review-nlfix-d7.R, gen(seed * 100 + 1),
# seeds 1:5, slopes +-0.6. Args: arm, part (fit | ref), group:
#   disc   cumulative, mu body x nlf(disc ~ exp(d1 * d2 * z))
#   sigma  gaussian, mu body x nlf(sigma ~ s0 + s1 * s2 * z) (log link)
#   mv     mvbf of two cumulative bodies b1 exp(b2 x) and c1 exp(c2 z)
# ref: the best of 8 random starts of every coefficient, set.seed(seed * 7 + 1).
a <- commandArgs(TRUE)
arm <- a[1]; part <- a[2]; grp <- a[3]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
mub <- list(exp = "b1 * exp(b2 * x)", sat = "b1 * (1 - exp(-b2 * x))",
            two = "b1 * exp(b2 * x) + b3 * z")
gen <- function(seed, slope) {
  set.seed(seed); n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  lat <- 1.5 * exp(slope * d$x) + 0.3 * d$z
  d$yo <- ordered(cut(lat * exp(0.3 * d$z) + rlogis(n),
                      c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  d$yg <- lat + rnorm(n, 0, 0.5 * exp(0.4 * d$z))
  d$y2 <- ordered(cut(1.2 * exp(-0.5 * d$z) + rlogis(n),
                      c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  d
}
pars_of <- function(s) setdiff(all.vars(str2lang(s)), c("x", "z"))
model <- function(bn) {
  if (grp == "disc") {
    p <- c(pars_of(mub[[bn]]), "d1", "d2")
    f <- bf(yo ~ 1) + nlf(as.formula(paste("mu ~", mub[[bn]]))) +
      nlf(disc ~ exp(d1 * d2 * z))
    for (q in p) f <- f + lf(as.formula(paste(q, "~ 1")))
    list(f = f, fam = cumulative(), p = p)
  } else if (grp == "sigma") {
    p <- c(pars_of(mub[[bn]]), "s0", "s1", "s2")
    f <- bf(yg ~ 1) + nlf(as.formula(paste("mu ~", mub[[bn]]))) +
      nlf(sigma ~ s0 + s1 * s2 * z)
    for (q in p) f <- f + lf(as.formula(paste(q, "~ 1")))
    list(f = f, fam = gaussian(), p = p)
  } else {
    f1 <- bf(yo ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE)
    f2 <- bf(y2 ~ c1 * exp(c2 * z), c1 + c2 ~ 1, nl = TRUE)
    list(f = mvbf(f1, f2, rescor = FALSE), fam = cumulative(),
         p = c("yo_b1", "yo_b2", "y2_c1", "y2_c2"))
  }
}
out <- file.path(FE_OUT, sprintf("nl2-%s-%s-%s.tsv", grp, part, arm))
cat("grp\tbody\tseed\tslope\tll\tconv\tnwarn\tsec\n", file = out)
bns <- if (grp == "mv") "mv" else names(mub)
for (seed in 1:5) for (slope in c(0.6, -0.6)) for (bn in bns) {
  d <- gen(seed * 100 + 1, slope)
  m <- model(bn)
  t0 <- proc.time()[[3]]
  if (part == "fit") {
    r <- suppressMessages(fe_fit(m$f, family = m$fam, data = d))
    ll <- if (inherits(r$fit, "error")) NA else as.numeric(logLik(r$fit))
    cv <- if (inherits(r$fit, "error")) NA else r$fit$opt$convergence
    nw <- length(r$warn)
  } else {
    set.seed(seed * 7 + 1)
    ll <- -Inf
    tp <- par_template(m$f, family = m$fam, data = d)
    for (i in 1:8) {
      st <- list()
      for (cp in c("beta", "betad")) {
        v <- tp[[cp]]
        if (length(v)) st[[cp]] <- setNames(rnorm(length(v), 0, 1), names(v))
      }
      r <- suppressMessages(fe_fit(m$f, family = m$fam, data = d, start = st))
      if (!inherits(r$fit, "error")) {
        v <- as.numeric(logLik(r$fit))
        if (is.finite(v) && v > ll) ll <- v
      }
    }
    cv <- NA; nw <- NA
  }
  cat(sprintf("%s\t%s\t%d\t%.1f\t%.6f\t%s\t%s\t%.2f\n", grp, bn, seed, slope, ll,
              cv, nw, proc.time()[[3]] - t0), file = out, append = TRUE)
}
cat("DONE", out, "\n")
