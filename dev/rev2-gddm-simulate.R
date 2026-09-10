# rev-gddm, second pass: the false-alarm rate of the new
# gddm_simulate() refusal.
#
# A simulator that refuses a legitimate call costs what a fit refusing
# one costs. The lane measured four legitimate shapes; this measures
# sixteen, plus five that must be refused, on the same footing as the
# density check's paired sweep: an OK arm that must stay silent and a
# BAD arm that must fire and must name the parameter.
#
# Seed 606. Arm from GDDM_LIB; default is this review's own install.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/rev-gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, " frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n\n")

set.seed(606)
ct <- gddm_control(t_max = 2)
n <- 100L
h2 <- rep(0:1, each = n / 2L)
h4 <- rep(1:4, each = n / 4L)

ok <- list()
add <- function(label, expr) ok[[label]] <<- substitute(expr)

add("scalars only",
    gddm_simulate(n, mu = 2, bs = 1.5, ndt = 0.2, control = ct))
add("length-n mu constant within two coh levels",
    gddm_simulate(n, mu = ifelse(h2 == 0, 1, 3), bs = 1.5, ndt = 0.2,
                  coh = h2, control = ct))
add("length-n mu constant within four coh levels",
    gddm_simulate(n, mu = c(1, 2, 3, 4)[h4], bs = 1.5, ndt = 0.2,
                  coh = h4, control = ct))
add("length-2 mu recycled, coh alternating with it",
    gddm_simulate(n, mu = c(1, 3), bs = 1.5, ndt = 0.2,
                  coh = rep(0:1, length.out = n), control = ct))
add("coh given as a character vector",
    gddm_simulate(n, mu = ifelse(h2 == 0, 1, 3), bs = 1.5, ndt = 0.2,
                  coh = ifelse(h2 == 0, "lo", "hi"), control = ct))
add("coh given as a factor",
    gddm_simulate(n, mu = ifelse(h2 == 0, 1, 3), bs = 1.5, ndt = 0.2,
                  coh = factor(h2), control = ct))
add("every parameter a vector, all constant within coh",
    gddm_simulate(n, mu = ifelse(h2 == 0, 1, 3),
                  bs = ifelse(h2 == 0, 1.5, 2), ndt = rep(0.2, n),
                  bias = rep(0.5, n), coh = h2, control = ct))
add("ndt a vector constant within coh",
    gddm_simulate(n, mu = 2, bs = 1.5, ndt = ifelse(h2 == 0, .2, .3),
                  coh = h2, control = ct))
add("a computed within-coh constant, ave()",
    gddm_simulate(n, mu = ave(rnorm(n), h2, FUN = mean) + 2, bs = 1.5,
                  ndt = 0.2, coh = h2, control = ct))
add("lapse uniform, lapse scalar",
    gddm_simulate(n, mu = 2, bs = 1.5, ndt = 0.2, lapse = 0.02,
                  control = gddm_control(t_max = 2)))
add("collapsing bound, tau constant within coh",
    gddm_simulate(n, mu = 2, bs = 1.5, ndt = 0.2,
                  tau = ifelse(h2 == 0, 1, 1.5), coh = h2,
                  bound = gddm_bound_exponential(), control = ct))
add("uniform start, sz constant within coh",
    gddm_simulate(n, mu = 2, bs = 1.5, ndt = 0.2,
                  sz = ifelse(h2 == 0, 0.05, 0.1), coh = h2,
                  start = gddm_start_uniform(), control = ct))
add("leaky drift, leak scalar",
    gddm_simulate(n, mu = 2, bs = 1.5, ndt = 0.2, leak = 0.5,
                  drift = list(gddm_drift_constant(),
                               gddm_drift_leak()), control = ct))
add("coherence drift, coh the covariate it reads",
    gddm_simulate(n, mu = 2, bs = 1.5, ndt = 0.2, alpha = 1,
                  coh = ifelse(h2 == 0, 0.2, 0.6),
                  drift = gddm_drift_coherence(), control = ct))
add("n = 1", gddm_simulate(1L, mu = 2, bs = 1.5, ndt = 0.2,
                           control = ct))
add("one coh level, one trial per parameter setting",
    gddm_simulate(4L, mu = c(1, 2, 3, 4), bs = 1.5, ndt = 0.2,
                  coh = 1:4, control = ct))
add("the help page's own example",
    gddm_simulate(20, mu = 1.5, bs = 2, ndt = 0.2,
                  control = gddm_control(t_max = 2, dt = 0.02,
                                         ny = 101)))

bad <- list()
addb <- function(label, expr, who) bad[[label]] <<-
  list(e = substitute(expr), who = who)
addb("mu varying inside one coh level",
     gddm_simulate(n, mu = c(rep(-2.5, n / 2L), rep(2.5, n / 2L)),
                   bs = 1.5, ndt = 0.2, coh = 0, control = ct),
     "mu")
addb("mu drawn per trial",
     gddm_simulate(n, mu = rnorm(n, 2), bs = 1.5, ndt = 0.2,
                   control = ct), "mu")
addb("ndt varying inside a coh level",
     gddm_simulate(n, mu = 2, bs = 1.5,
                   ndt = rep(c(0.2, 0.3), length.out = n), coh = h2,
                   control = ct), "ndt")
addb("two parameters varying",
     gddm_simulate(n, mu = rnorm(n, 2), bs = runif(n, 1, 2),
                   ndt = 0.2, control = ct), "bs")
addb("bias varying inside a coh level",
     gddm_simulate(n, mu = 2, bs = 1.5, ndt = 0.2,
                   bias = rep(c(0.4, 0.6), length.out = n), coh = h2,
                   control = ct), "bias")

cat("== OK arm: legitimate calls, which must NOT be refused\n")
nfa <- 0L
for (nm in names(ok)) {
  r <- tryCatch({ z <- eval(ok[[nm]]); sprintf("ok, %d rows", nrow(z)) },
                error = function(e) paste("REFUSED:",
                                          conditionMessage(e)))
  if (grepl("^REFUSED", r)) nfa <- nfa + 1L
  cat(sprintf("%s %-46s %s\n", if (grepl("^REFUSED", r)) ">>>" else "  ",
              nm, substr(r, 1, 90)))
}
cat("\n  legitimate shapes:", length(ok), "  FALSE ALARMS:", nfa, "\n")

cat("\n== BAD arm: calls that must be refused, and must name the\n")
cat("   parameter\n")
nmiss <- 0L; nwrong <- 0L
for (nm in names(bad)) {
  b <- bad[[nm]]
  r <- tryCatch({ eval(b$e); "" },
                error = function(e) conditionMessage(e))
  if (!nzchar(r)) {
    nmiss <- nmiss + 1L
    cat(sprintf(">>> %-40s MISSED\n", nm))
  } else if (!grepl(paste0("`", b$who, "`"), r, fixed = TRUE)) {
    nwrong <- nwrong + 1L
    cat(sprintf(">>> %-40s named wrong: %s\n", nm, substr(r, 1, 70)))
  } else {
    cat(sprintf("    %-40s refused, names `%s`\n", nm, b$who))
  }
}
cat("\n  bad shapes:", length(bad), "  missed:", nmiss,
    "  unnamed:", nwrong, "\n")
