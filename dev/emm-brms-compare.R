# emmeans on brms and frmtmb, same data, same arguments (lane emm).
#
#   Rscript dev/emm-brms-compare.R > dev/emm-brms-compare-log.txt 2>&1
#
# Part 1 needs no sampler: brms's recover_data() on an `empty = TRUE`
# fit shows which variables brms puts in the reference grid and which
# argument combinations it refuses. Part 2 samples two small models
# (flat priors, 4 x 1000 draws) and prints brms's marginal means and
# posterior SDs beside frmtmb's estimates and standard errors. The two
# are different estimators, so the comparison is of meaning and layout
# (which predictor, which scale, which levels), with agreement to about
# the Monte Carlo error expected, not bit-for-bit.
#
# Part 2 runs only with EMM_BRMS_SAMPLE=true. On the lane-emm machine
# (2026-09-25) rstan could not compile any model: StanHeaders includes
# tbb/tbb_stddef.h, which the conda RcppParallel build does not ship.
# The log therefore holds part 1 only.
lib <- Sys.getenv("FRMTMB_LIB", "/opt/rlib/lane-emm")
.libPaths(c(lib, "/opt/rlib/base", "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages({
  library(brms)
  library(emmeans)
})
set.seed(1)
n <- 200
d <- data.frame(f = factor(sample(c("A", "B", "C"), n, TRUE)),
                x = runif(n), g = factor(sample(1:6, n, TRUE)))
d$y1 <- 1 + as.numeric(d$f) * 0.5 + 2 * d$x + rnorm(n)
d$y2 <- 0.5 * as.numeric(d$f) + rnorm(n)

show <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.data.frame(r)) {
    cat(sprintf("%-28s grid variables: %s\n", label,
                paste(attr(r, "predictors"), collapse = " ")))
  } else {
    cat(sprintf("%-28s %s\n", label, gsub("\n", " ", r)))
  }
}

cat("== 1. brms's reference-grid variables and refusals\n")
enl <- brm(bf(y1 ~ a + b * x, a ~ f + (1 | g), b ~ 1, nl = TRUE),
           data = d, empty = TRUE)
show("nl, no nlpar", recover_data(enl))
show("nl, nlpar = 'a'", recover_data(enl, nlpar = "a"))
show("nl, nlpar a, re NULL", recover_data(enl, nlpar = "a",
                                          re_formula = NULL))
show("nl, epred", recover_data(enl, epred = TRUE))
show("nl, dpar + nlpar", recover_data(enl, dpar = "sigma", nlpar = "a"))
show("nl, unknown nlpar", recover_data(enl, nlpar = "zz"))
emv <- brm(bf(y1 ~ f) + bf(y2 ~ x) + set_rescor(FALSE), data = d,
           empty = TRUE)
show("mv, no resp", recover_data(emv))
show("mv, resp = 'y1'", recover_data(emv, resp = "y1"))
show("mv, unknown resp", recover_data(emv, resp = "y3"))
d$cnt <- rpois(n, 3)
emx <- brm(bf(y1 ~ f) + bf(cnt ~ f, family = poisson()) +
             set_rescor(FALSE), data = d, empty = TRUE)
show("mv, mixed families", recover_data(emx))
show("mv, mixed families, epred", recover_data(emx, epred = TRUE))

cat("\n== 1b. frmtmb, the same calls\n")
suppressMessages(library(frmtmb))
fnl1 <- frm(frmtmb::bf(y1 ~ a + b * x, a ~ f + (1 | g), b ~ 1, nl = TRUE),
            data = d)
show("nl, no nlpar", recover_data(fnl1))
show("nl, nlpar = 'a'", recover_data(fnl1, nlpar = "a"))
show("nl, nlpar a, re NULL", recover_data(fnl1, nlpar = "a",
                                          re_formula = NULL))
show("nl, epred", recover_data(fnl1, epred = TRUE))
show("nl, dpar + nlpar", recover_data(fnl1, dpar = "sigma", nlpar = "a"))
show("nl, unknown nlpar", recover_data(fnl1, nlpar = "zz"))
fmv1 <- frm(mvbf(frmtmb::bf(y1 ~ f), frmtmb::bf(y2 ~ x)), data = d)
show("mv, no resp", recover_data(fmv1))
show("mv, resp = 'y1'", recover_data(fmv1, resp = "y1"))
show("mv, unknown resp", recover_data(fmv1, resp = "y3"))
fmx1 <- frm(mvbf(frmtmb::bf(y1 ~ f), frmtmb::bf(cnt ~ f) + poisson()),
            data = d)
show("mv, mixed families", recover_data(fmx1))
show("mv, mixed families, epred", recover_data(fmx1, epred = TRUE))
cat("rep.meas levels, frmtmb:",
    levels(ref_grid(fmv1)@grid$rep.meas), "\n")

if (!identical(Sys.getenv("EMM_BRMS_SAMPLE"), "true")) {
  cat("\npart 2 skipped: set EMM_BRMS_SAMPLE=true to sample\n")
  quit(save = "no")
}
cat("\n== 2. sampled models against frmtmb\n")
suppressMessages(library(frmtmb))
side <- function(label, eb, ef) {
  sb <- as.data.frame(summary(eb, point.est = mean))
  sf <- as.data.frame(summary(ef))
  pb <- as.matrix(as.mcmc(eb))
  cat(label, "\n")
  print(data.frame(brms_mean = round(sb$emmean, 4),
                   brms_sd = round(apply(pb, 2, sd), 4),
                   frmtmb_est = round(sf$emmean, 4),
                   frmtmb_se = round(sf$SE, 4),
                   row.names = do.call(paste, sf[seq_len(
                     which(names(sf) == "emmean") - 1L)])))
}
# frmtmb masks bf(), prior() and set_rescor(), so part 2 names brms's
bnl <- brm(brms::bf(y1 ~ a + b * x, a ~ f, b ~ 1, nl = TRUE), data = d,
           prior = c(brms::prior(normal(0, 100), nlpar = "a"),
                     brms::prior(normal(0, 100), nlpar = "b")),
           chains = 4, iter = 2000, seed = 1, refresh = 0)
fnl <- frm(bf(y1 ~ a + b * x, a ~ f, b ~ 1, nl = TRUE), data = d)
side("nlpar = 'a'", emmeans(bnl, "f", nlpar = "a"),
     emmeans(fnl, "f", nlpar = "a"))
side("epred = TRUE", emmeans(bnl, "f", epred = TRUE),
     emmeans(fnl, "f", epred = TRUE))
side("pairs, epred = TRUE", pairs(emmeans(bnl, "f", epred = TRUE)),
     pairs(emmeans(fnl, "f", epred = TRUE)))
cat("brms, no nlpar, over x only:\n")
print(summary(emmeans(bnl, "x")))
cat("frmtmb, no nlpar, over x (f averaged):\n")
print(summary(emmeans(fnl, "x")))

bmv <- brm(brms::bf(y1 ~ f) + brms::bf(y2 ~ f) + brms::set_rescor(FALSE),
           data = d,
           chains = 4, iter = 2000, seed = 1, refresh = 0)
fmv <- frm(mvbf(bf(y1 ~ f), bf(y2 ~ f)), data = d)
eb <- emmeans(bmv, ~ f | rep.meas)
ef <- emmeans(fmv, ~ f | rep.meas)
cat("rep.meas levels: brms", levels(eb@grid$rep.meas), "| frmtmb",
    levels(ef@grid$rep.meas), "\n")
side("no resp, f | rep.meas", eb, ef)
side("resp = 'y2'", emmeans(bmv, "f", resp = "y2"),
     emmeans(fmv, "f", resp = "y2"))
