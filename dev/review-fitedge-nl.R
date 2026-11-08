# Review of lane fitedge, item 2. Args: arm (base | lane), part (fit | ref),
# group:
#   one   24 single-body models (6 families x 4 bodies); every value written in
#         hex, for a bitwise comparison of the arms
#   tri   skew_normal with three bodies: mu, sigma and alpha
#   mv4   mvbf of a cumulative (mu and disc bodies) and a gaussian (mu and
#         sigma bodies) response: four bodies
#   mvs   mvbf of two gaussian responses, mu and sigma bodies in each
# Data: gen(seed * 100 + 7, slope), seeds 1:4, slopes +-0.6.
# ref: the best of 8 random starts of every coefficient, set.seed(seed * 13 + 5).
a <- commandArgs(TRUE)
arm <- a[1]; part <- a[2]; grp <- a[3]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
RV_OUT <- file.path(FE_WT, "dev/review-fitedge-out")
gen <- function(seed, slope) {
  set.seed(seed); n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3))
  lat <- 1.5 * exp(slope * d$x) + 0.3 * d$z
  d$yo <- ordered(cut(lat * exp(0.3 * d$z) + rlogis(n),
                      c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))
  d$yg <- lat + rnorm(n, 0, 0.5 * exp(0.4 * d$z))
  d$ys <- lat + 0.5 * exp(0.3 * d$z) * (rexp(n) - 1) * (1 + 0.5 * (d$w > 0))
  d$yk <- rpois(n, exp(0.3 + 0.5 * exp(slope * d$x) / 2))
  d$yb <- rbinom(n, 1, plogis(-1 + exp(slope * d$x)))
  d$ya <- rgamma(n, 4, 4 / exp(0.2 + 0.4 * exp(slope * d$x) / 2))
  d$yt <- lat + 0.5 * rt(n, 4)
  d$yg2 <- 0.8 * exp(-slope * d$z) + rnorm(n, 0, 0.4 * exp(-0.3 * d$x))
  d
}
pars_of <- function(s) setdiff(all.vars(str2lang(s)), c("x", "z", "w"))
add_lf <- function(f, p) {
  for (q in p) f <- f + lf(as.formula(paste(q, "~ 1")))
  f
}
one_body <- list(exp = "b1 * exp(b2 * x)", sat = "b1 * (1 - exp(-b2 * x))",
                 two = "b1 * exp(b2 * x) + b3 * z", lgt = "b1 / (1 + exp(-b2 * x))")
one_fam <- list(gaussian = list("yg", gaussian()), poisson = list("yk", poisson()),
                bernoulli = list("yb", bernoulli()),
                cumulative = list("yo", cumulative()),
                gamma = list("ya", Gamma(link = "log")),
                student = list("yt", student()))
models <- function(grp) {
  out <- list()
  if (grp == "one") {
    for (fm in names(one_fam)) for (bn in names(one_body)) {
      rsp <- one_fam[[fm]][[1]]
      p <- pars_of(one_body[[bn]])
      f <- add_lf(bf(as.formula(paste(rsp, "~ 1"))) +
                    nlf(as.formula(paste("mu ~", one_body[[bn]]))), p)
      out[[paste(fm, bn, sep = "_")]] <- list(f = f, fam = one_fam[[fm]][[2]])
    }
  } else if (grp == "tri") {
    for (bn in c("exp", "sat")) {
      p <- c(pars_of(one_body[[bn]]), "s0", "s1", "s2", "a1", "a2")
      f <- add_lf(bf(ys ~ 1) + nlf(as.formula(paste("mu ~", one_body[[bn]]))) +
                    nlf(sigma ~ s0 + s1 * s2 * z) + nlf(alpha ~ a1 * a2 * w), p)
      out[[bn]] <- list(f = f, fam = skew_normal())
    }
  } else if (grp == "mv4") {
    f1 <- add_lf(bf(yo ~ 1, family = cumulative()) + nlf(mu ~ b1 * exp(b2 * x)) +
                   nlf(disc ~ exp(d1 * d2 * z)), c("b1", "b2", "d1", "d2"))
    f2 <- add_lf(bf(yg ~ 1, family = gaussian()) + nlf(mu ~ c1 * exp(c2 * x)) +
                   nlf(sigma ~ s0 + s1 * s2 * z), c("c1", "c2", "s0", "s1", "s2"))
    out[["mv4"]] <- list(f = mvbf(f1, f2, rescor = FALSE), fam = NULL)
  } else if (grp == "mvs") {
    f1 <- add_lf(bf(yg ~ 1) + nlf(mu ~ b1 * exp(b2 * x)) +
                   nlf(sigma ~ s0 + s1 * s2 * z), c("b1", "b2", "s0", "s1", "s2"))
    f2 <- add_lf(bf(yg2 ~ 1) + nlf(mu ~ c1 * exp(c2 * z)) +
                   nlf(sigma ~ t0 + t1 * t2 * x), c("c1", "c2", "t0", "t1", "t2"))
    out[["mvs"]] <- list(f = mvbf(f1, f2, rescor = FALSE), fam = gaussian())
  }
  out
}
fitm <- function(m, d, ...) {
  if (is.null(m$fam)) {
    suppressMessages(fe_fit(m$f, data = d, ...))
  } else {
    suppressMessages(fe_fit(m$f, family = m$fam, data = d, ...))
  }
}
hex <- function(v) paste(sprintf("%a", as.numeric(v)), collapse = ",")
out <- file.path(RV_OUT, sprintf("nl-%s-%s-%s.tsv", grp, part, arm))
cat("grp\tmodel\tseed\tslope\tll\tconv\tnwarn\twarn\tsec\tparhex\tsehex\n",
    file = out)
ms <- models(grp)
seeds <- if (grp == "one") 1L else 1:4
slopes <- if (grp == "one") 0.6 else c(0.6, -0.6)
for (seed in seeds) for (slope in slopes) for (mn in names(ms)) {
  sl <- if (grp == "one" && match(mn, names(ms)) %% 2L == 0L) -0.6 else slope
  d <- gen(seed * 100 + 7, sl)
  m <- ms[[mn]]
  t0 <- proc.time()[[3]]
  ph <- sh <- wtxt <- ""
  if (part == "fit") {
    r <- fitm(m, d)
    ok <- !inherits(r$fit, "error")
    ll <- if (ok) as.numeric(logLik(r$fit)) else NA
    cv <- if (ok) r$fit$opt$convergence else NA
    nw <- length(r$warn)
    wtxt <- if (ok) paste(sort(substr(gsub("[\t\n]", " ", r$warn), 1, 90)),
                          collapse = " || ") else
      paste("ERROR", substr(conditionMessage(r$fit), 1, 90))
    if (ok) {
      ph <- hex(r$fit$opt$par)
      sh <- hex(tryCatch(suppressWarnings(sqrt(diag(as.matrix(
        vcov(r$fit, full = TRUE))))), error = function(e) NA))
    }
  } else {
    set.seed(seed * 13 + 5)
    ll <- -Inf
    tp <- if (is.null(m$fam)) par_template(m$f, data = d) else
      par_template(m$f, family = m$fam, data = d)
    for (i in 1:8) {
      st <- list()
      for (cp in c("beta", "betad")) {
        v <- tp[[cp]]
        if (length(v)) st[[cp]] <- setNames(rnorm(length(v), 0, 1), names(v))
      }
      r <- fitm(m, d, start = st)
      if (!inherits(r$fit, "error")) {
        v <- as.numeric(logLik(r$fit))
        if (is.finite(v) && v > ll) ll <- v
      }
    }
    cv <- NA; nw <- NA
  }
  cat(sprintf("%s\t%s\t%d\t%.1f\t%.8f\t%s\t%s\t%s\t%.2f\t%s\t%s\n", grp, mn,
              seed, sl, ll, cv, nw, wtxt, proc.time()[[3]] - t0, ph, sh),
      file = out, append = TRUE)
}
cat("DONE", out, "\n")
