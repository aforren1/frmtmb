# Does predict()'s interval cover at the nominal rate?
#
# The question is a PREDICTIVE one, so the measurement is out of sample:
# fit on n rows, form the interval at m FRESH design points, and draw
# the responses at those points from the same truth. Covering the
# training rows would ask a different question, because the fit has
# already chased them.
#
# The arms, so that the construction is justified rather than assumed:
#   joint    parameters drawn from N(theta_hat, vcov(full = TRUE))
#   plugin   simulated at the estimates alone (param_uncertainty = FALSE)
#   wald     the Wald interval around the FITTED value, which carries no
#            observation noise at all: the control that MUST under-cover
#            badly, so that a coverage near 0.95 in the other arms is
#            not an artifact of the harness
#   exact    (gauss only) the textbook prediction interval from lm(),
#            t_{n-p} sigma sqrt(1 + h): a positive control that MUST
#            cover at the nominal rate, so that the harness is shown to
#            be able to report 0.95 as well as to report a shortfall
#   condvar  (mixed only) the joint arm widened by the analytic
#            conditional variance of the group mode,
#            sigma^2 tau^2 / (sigma^2 + n_g tau^2). This is not a
#            proposed construction: it is the TEST of the attribution
#            the documentation makes, that the mixed shortfall is the
#            missing Var(b | y). If adding exactly that term closes the
#            gap, the attribution is measured; if it does not, the
#            attribution is wrong.
#
#   Rscript dev/shapes-coverage.R <design> <nrep>
#     > dev/shapes-log/cov-<design>.txt
# designs: gauss, pois, mixed, newlevel
#
# The `newlevel` design is the one BLOCKER 2 of
# dev/reviews/20260918-shapes.md was about: every fresh point sits at a
# grouping level the fit never saw, so the correct predictive spread is
# sqrt(sigma^2 + tau^2) and not sqrt(sigma^2). The truth draws a FRESH
# group effect per new level, which is what "a level the fit has not
# seen" means.
#
# Power. This is a ONE-sample count against a rate known exactly, so
# `power.prop.test` is the wrong instrument: 202 replicates is where
# 0.80 power is first SUSTAINED for a one-sample count at 0.95
# (dev/lane-rules.md). 220 replicates of 10 new points each is 2200
# predictions per arm. The predictions inside one replicate share a
# fit, so the pooled binomial interval is too narrow; the
# replicate-level mean and its standard error are reported beside it
# and are the honest ones.

.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))

a <- commandArgs(trailingOnly = TRUE)
design <- if (length(a)) a[1] else "gauss"
nrep <- if (length(a) > 1) as.integer(a[2]) else 220L
n <- 60L
m <- 10L
ndraws <- 1000L
seed0 <- 20260917L
TAU <- 0.7
SIGMA <- 1

gen <- function(design, n, seed) {
  set.seed(seed)
  if (design == "gauss") {
    d <- data.frame(x = rnorm(n))
    d$y <- rnorm(n, 1 + 0.5 * d$x, 1)
    list(data = d, form = bf(y ~ x), fam = gaussian(),
         newx = function(k) data.frame(x = rnorm(k)),
         draw = function(nd) rnorm(nrow(nd), 1 + 0.5 * nd$x, 1))
  } else if (design == "pois") {
    d <- data.frame(x = rnorm(n))
    d$y <- rpois(n, exp(0.4 + 0.5 * d$x))
    list(data = d, form = bf(y ~ x), fam = poisson(),
         newx = function(k) data.frame(x = rnorm(k)),
         draw = function(nd) rpois(nrow(nd), exp(0.4 + 0.5 * nd$x)))
  } else if (design == "mixed") {
    ng <- 12L
    d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng),
                                                 length.out = n)))
    u <- rnorm(ng, 0, TAU)
    d$y <- rnorm(n, 1 + 0.5 * d$x + u[d$g], SIGMA)
    list(data = d, form = bf(y ~ x + (1 | g)), fam = gaussian(),
         u = u, ng = ng, group = d$g,
         newx = function(k) {
           data.frame(x = rnorm(k),
                      g = factor(sample(seq_len(ng), k, TRUE),
                                 levels = levels(d$g)))
         },
         draw = function(nd) {
           rnorm(nrow(nd), 1 + 0.5 * nd$x + u[as.integer(nd$g)], SIGMA)
         })
  } else if (design == "newlevel") {
    ng <- 12L
    d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(ng),
                                                 length.out = n)))
    u <- rnorm(ng, 0, TAU)
    d$y <- rnorm(n, 1 + 0.5 * d$x + u[d$g], SIGMA)
    list(data = d, form = bf(y ~ x + (1 | g)), fam = gaussian(),
         u = u, ng = ng, new_levels = TRUE,
         newx = function(k) {
           # levels the fit has never seen, one per fresh point
           data.frame(x = rnorm(k), g = factor(paste0("new", seq_len(k))))
         },
         draw = function(nd) {
           # a fresh group effect per unseen level, which is what the
           # model says an unseen level is
           unew <- rnorm(nrow(nd), 0, TAU)
           rnorm(nrow(nd), 1 + 0.5 * nd$x + unew, SIGMA)
         })
  } else {
    stop("unknown design")
  }
}

arms <- switch(design,
               gauss = c("joint", "plugin", "wald", "exact"),
               mixed = c("joint", "plugin", "wald", "condvar"),
               c("joint", "plugin", "wald"))
hits <- stats::setNames(rep(list(integer(0)), length(arms)), arms)
widths <- stats::setNames(rep(list(numeric(0)), length(arms)), arms)
tot <- 0L
fail <- 0L
anl <- isTRUE(gen(design, n, seed0)$new_levels)
z <- stats::qnorm(0.975)

for (r in seq_len(nrep)) {
  seed <- seed0 + r
  g <- gen(design, n, seed)
  fit <- tryCatch(suppressWarnings(frm(g$form + g$fam, data = g$data)),
                  error = function(e) NULL)
  if (is.null(fit)) { fail <- fail + 1L; next }
  nd <- g$newx(m)
  ynew <- g$draw(nd)
  set.seed(seed * 7L)
  pj <- tryCatch(suppressWarnings(
    predict(fit, newdata = nd, ndraws = ndraws,
            allow_new_levels = anl)), error = function(e) NULL)
  set.seed(seed * 7L + 1L)
  pp <- tryCatch(suppressWarnings(
    predict(fit, newdata = nd, ndraws = ndraws, allow_new_levels = anl,
            param_uncertainty = FALSE)), error = function(e) NULL)
  pw <- tryCatch(suppressWarnings(
    fitted(fit, newdata = nd, allow_new_levels = anl)),
    error = function(e) NULL)
  if (is.null(pj) || is.null(pp) || is.null(pw)) {
    fail <- fail + 1L
    next
  }
  got <- list(joint = pj, plugin = pp, wald = pw)
  if ("exact" %in% arms) {
    lmf <- stats::lm(y ~ x, data = g$data)
    pe <- stats::predict(lmf, newdata = nd, interval = "prediction")
    got$exact <- cbind(pe[, "fit"], NA, pe[, "lwr"], pe[, "upr"])
  }
  if ("condvar" %in% arms) {
    # Var(b | y) for a balanced one-way intercept block, at the
    # ESTIMATES: sigma^2 tau^2 / (sigma^2 + n_g tau^2)
    s2 <- sigma(fit)^2
    t2 <- VarCorr(fit)$g$sd[1L, "Estimate"]^2
    ng_obs <- as.numeric(table(g$group)[as.character(nd$g)])
    v <- s2 * t2 / (s2 + ng_obs * t2)
    sd_j <- (pj[, 4L] - pj[, 3L]) / (2 * z)
    wide <- sqrt(sd_j^2 + v)
    got$condvar <- cbind(pj[, 1L], NA, pj[, 1L] - z * wide,
                         pj[, 1L] + z * wide)
  }
  for (nm in arms) {
    p <- got[[nm]]
    hits[[nm]] <- c(hits[[nm]], sum(ynew >= p[, 3L] & ynew <= p[, 4L]))
    widths[[nm]] <- c(widths[[nm]], mean(p[, 4L] - p[, 3L]))
  }
  tot <- tot + m
}

cat("design:", design, " replicates:", length(hits$joint),
    " failed fits:", fail, "\n")
cat("n:", n, " new points per replicate:", m, " ndraws:", ndraws,
    " seed0:", seed0, " allow_new_levels:", anl, "\n")
cat("script: dev/shapes-coverage.R\n\n")
for (nm in arms) {
  k <- sum(hits[[nm]])
  ci <- stats::binom.test(k, tot, 0.95)$conf.int
  # the replicate-level standard error, which the pooled binomial one
  # understates because the m predictions of one replicate share a fit
  pr <- hits[[nm]] / m
  se_r <- stats::sd(pr) / sqrt(length(pr))
  cat(sprintf("%-8s %5d of %5d = %.4f  (%.4f, %.4f) binomial; ",
              nm, k, tot, k / tot, ci[1], ci[2]))
  # the MEDIAN width, not the mean: a parametric-bootstrap draw of a
  # log-scale parameter can overflow to Inf on one replicate in a
  # thousand, which makes a mean width uninformative. The count of
  # such replicates is printed rather than dropped in silence.
  w <- widths[[nm]]
  cat(sprintf("replicate mean %.4f se %.4f  median width %.4f",
              mean(pr), se_r, stats::median(w[is.finite(w)])))
  cat(sprintf("  non-finite widths %d of %d\n",
              sum(!is.finite(w)), length(w)))
}
