# rev-gddm: constructions where the guarded thing is ABSENT, at frame
# assembly only (dry_run = "frame"), so nothing is taped.
#
# Two questions:
#   1. can a construction make the refusal fall silent on a model whose
#      dpar really does vary inside a condition?
#   2. can a construction make it fire on a model that is correct?
#
# Seed 909. Arm from GDDM_LIB; default is this review's own install of
# the worktree, C:/Users/adf44/source/r/rev-gddm-lib.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/rev-gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, " frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n\n")

ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
set.seed(909)
n <- 120L
d0 <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                    control = gddm_control(t_max = 2))

# helper: run one frame assembly and report the refusal, if any
run <- function(data, spec, fam = gddm(control = ctl)) {
  tryCatch({
    frm(spec, family = fam, data = data, dry_run = "frame")
    ""
  }, error = function(e) conditionMessage(e))
}
say <- function(label, msg, want) {
  got <- if (nzchar(msg)) "REFUSED" else "accepted"
  flag <- if (got == want) "   " else ">>>"
  cat(sprintf("%s %-52s %s\n", flag, label, got))
  if (nzchar(msg)) cat("      ", substr(gsub("\\s+", " ", msg), 1, 96),
                       "\n", sep = "")
  invisible(got == want)
}
mu_x <- bf(rt | vint(upper, cond) ~ x, bs ~ 1, ndt ~ 1, bias = 0.5)

cat("== A. the guard PRESENT: the base case both arms agree on\n")
d <- d0
d$cond <- rep(1:2, length.out = n)      # 2 conditions, 60 rows each
d$x <- rnorm(n)                          # varies inside both
say("mu ~ x, x varying inside a condition", run(d, mu_x), "REFUSED")
d$xc <- ifelse(d$cond == 1L, 0.3, 0.7)   # constant inside
say("mu ~ xc, xc constant inside a condition",
    run(d, bf(rt | vint(upper, cond) ~ xc, bs ~ 1, ndt ~ 1, bias = 0.5)),
    "accepted")

cat("\n== B. non-finite and degenerate columns\n")
d <- d0; d$cond <- rep(1:2, length.out = n); d$x <- rnorm(n)
b <- d; b$x[3] <- Inf
say("x varies, one Inf elsewhere in the column", run(b, mu_x), "REFUSED")
b <- d; b$x[3] <- NaN
say("x varies, one NaN elsewhere in the column", run(b, mu_x), "REFUSED")
b <- d; b$x <- rep(c(0, 0), length.out = n); b$x[4] <- 1e-300
say("all-zero column, one 1e-300 inside a condition",
    run(b, mu_x), "REFUSED")
b <- d; b$x <- rep(Inf, n)
say("column of all Inf (no variation)", run(b, mu_x), "accepted")
b <- d; b$x <- rep(Inf, n); b$x[2] <- 1
say("column of Inf with one finite entry inside", run(b, mu_x),
    "REFUSED")
b <- d; b$x <- rep(NA_real_, n)
say("column of all NA", run(b, mu_x), "accepted")
b <- d; b$x <- rnorm(n); b$x[2] <- NA
say("x varies and carries one NA (rows dropped by na.omit)",
    run(b, mu_x), "REFUSED")

cat("\n== C. the condition index itself\n")
b <- d; b$cond[5] <- NA_integer_
say("index with one NA", run(b, mu_x), "REFUSED")
b <- d; b$cond <- seq_len(n)
say("one row per condition, x varying between rows", run(b, mu_x),
    "accepted")
b <- d; b$cond <- c(1L, rep(2L, n - 1L))
say("one condition of 1 row, one of 119, x varying",
    run(b, mu_x), "REFUSED")
b <- d; b$cond <- rep(1:2, length.out = n); b$x <- rnorm(n)
b$x[b$cond == 2L] <- 0.7   # varies in condition 1 only
say("x varies inside condition 1 only", run(b, mu_x), "REFUSED")

cat("\n== D. factors\n")
b <- d; b$g <- factor(rep(c("a", "b"), length.out = n),
                      levels = c("a", "b", "unused"))
b$g[] <- "a"; b$g <- factor(b$g, levels = c("a", "b", "unused"))
say("factor constant, two unused levels",
    run(b, bf(rt | vint(upper, cond) ~ g, bs ~ 1, ndt ~ 1, bias = 0.5)),
    "accepted")
b <- d; b$g <- factor(ifelse(b$cond == 1L, "a", "b"),
                      levels = c("a", "b", "unused"))
say("factor constant within condition, one unused level",
    run(b, bf(rt | vint(upper, cond) ~ g, bs ~ 1, ndt ~ 1, bias = 0.5)),
    "accepted")
b <- d; b$g <- factor(sample(c("a", "b"), n, replace = TRUE))
say("factor varying within condition",
    run(b, bf(rt | vint(upper, cond) ~ g, bs ~ 1, ndt ~ 1, bias = 0.5)),
    "REFUSED")

cat("\n== E. the tolerance band, on covariates a field writes\n")
# The tolerance is 1e-8 times the column's own largest FINITE entry, so
# a column with a large additive offset or a long tail carries a
# tolerance that is large in the units the user reads.
band <- function(label, xv, want) {
  b <- d; b$x <- xv
  gi <- b$cond
  tol <- 1e-8 * max(abs(xv[is.finite(xv)]))
  spread <- max(tapply(xv, gi, function(z) max(z) - min(z)))
  cat(sprintf("   %-44s tol %.4g  within-condition spread %.4g\n",
              label, tol, spread))
  say(paste0("   -> ", label), run(b, mu_x), want)
}
# epoch milliseconds, as jsPsych and PsychoPy record a trial onset
t0 <- 1757462400000              # 2025-09-10 00:00:00 UTC, in ms
band("trial onset, epoch ms, 5 s inside a condition",
     t0 + rep(c(0, 5000), length.out = n) +
       rep(seq(0, 59) * 1000, each = 2), "REFUSED")
band("trial onset, epoch ms, 2 s inside a condition",
     t0 + rep(c(0, 2000), length.out = n), "REFUSED")
band("trial onset, epoch s, 2 s inside a condition",
     t0 / 1000 + rep(c(0, 2), length.out = n), "REFUSED")
# a delay-discounting covariate: seconds, from one second to ten years
band("delay in seconds, 1 s to 10 years, 1 s inside",
     ifelse(d$cond == 1L, rep(c(1, 2), length.out = n), 3.15e8),
     "REFUSED")

cat("\n== F. what the band accepts, priced in the user's own units\n")
for (tag in c("epoch ms", "epoch s", "delay s")) {
  xv <- switch(tag,
    "epoch ms" = t0 + rep(c(0, 2000), length.out = n),
    "epoch s"  = t0 / 1000 + rep(c(0, 2), length.out = n),
    "delay s"  = ifelse(d$cond == 1L, rep(c(1, 2), length.out = n),
                        3.15e8))
  mx <- max(abs(xv))
  cat(sprintf("   %-10s column max %.4g  tolerance %.4g %s\n",
              tag, mx, 1e-8 * mx,
              switch(tag, "epoch ms" = "ms", "epoch s" = "s",
                     "delay s" = "s")))
}
