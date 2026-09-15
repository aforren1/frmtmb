# Lane eamhier, item 2.1 of dev/extension-gaps-plan.md.
#
# One place for the design, the truths and the per-replicate record, so
# that every arm of the study draws from the same code and a number in
# dev/eamhier-findings.md can be re-derived from its seed alone.
#
# Sourced by eamhier-rep.R (one replicate, one process). Never run on
# its own.

eamhier_libs <- function() {
  .libPaths(c("C:/Users/adf44/source/r/eamhier-lib",
              "C:/Users/adf44/source/r/rellib-r3",
              "C:/Users/adf44/source/r/pinlib",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
}

# The tier's own truths (tests/testthat/test-scale.R, `eam_truth`), so
# that replicate 1 of arm A at seed 20260908 reproduces the recorded
# scale row rather than resembling it.
eamhier_truth <- list(mu0 = 0.4, mu_cond = 0.9, bs = 1.4, ndt = 0.25,
                      sd_mu = 0.35, sd_lbs = 0.20, sd_lndt = 0.12,
                      sv = 0.4)

# ------------------------------------------------------------ the data

# Arms A and C: the plan's realistic scale, 30 subjects x 400 trials, two
# conditions, random effects on the drift, the boundary and the
# non-decision time. Byte-for-byte the construction in
# tests/testthat/test-scale.R, with the seed as an argument.
#
# Arm B: the falsification design of dev/ndt-findings.md, "But the
# recovered spread is the FLOORS, not the random effect". `ndt` is
# CONSTANT across subjects and the log boundary spread is 0.25, so the
# per-subject floors vary and the truth does not. No condition effect
# and no drift random effect, because that is the design the recorded
# 6.4 ms and 23.6 ms came from.
eamhier_data <- function(seed, arm, ns = 30L, nt = 400L) {
  tr <- eamhier_truth
  set.seed(seed)
  if (identical(arm, "B")) {
    u_nd <- stats::rnorm(ns, 0, 0)
    u_bs <- stats::rnorm(ns, 0, 0.25)
    s <- rep(seq_len(ns), each = nt)
    t0 <- 0.25 * exp(u_nd)
    d <- frmtmb.eam::ddm_simulate(ns * nt, mu = 1.1,
                                  bs = 1.4 * exp(u_bs),
                                  ndt = t0[s], bias = 0.5)
    d$s <- factor(s)
    d$cond <- factor(rep(rep(0:1, each = nt / 2L), times = ns),
                     labels = c("a", "b"))
    attr(d, "ndt_subject") <- t0
    attr(d, "mu_subject") <- rep(1.1, ns)
    attr(d, "bs_subject") <- 1.4 * exp(u_bs)
    return(d)
  }
  sv <- if (identical(arm, "C")) tr$sv else 0
  u_mu <- stats::rnorm(ns, 0, tr$sd_mu)
  u_bs <- stats::rnorm(ns, 0, tr$sd_lbs)
  u_nd <- stats::rnorm(ns, 0, tr$sd_lndt)
  s <- rep(seq_len(ns), each = nt)
  cond <- rep(rep(0:1, each = nt / 2L), times = ns)
  d <- frmtmb.eam::ddm_simulate(ns * nt,
                                mu = tr$mu0 + tr$mu_cond * cond +
                                  u_mu[s],
                                bs = tr$bs * exp(u_bs[s]),
                                ndt = tr$ndt * exp(u_nd[s]),
                                bias = 0.5, sv = sv)
  d$s <- factor(s)
  d$cond <- factor(cond, labels = c("a", "b"))
  attr(d, "ndt_subject") <- tr$ndt * exp(u_nd)
  attr(d, "mu_subject") <- tr$mu0 + u_mu
  attr(d, "bs_subject") <- tr$bs * exp(u_bs)
  d
}

# `bound` is "pg" for the per-subject bound item 1.0a added and "gl" for
# 0.6.0's single global one. Arm B carries no condition effect and no
# drift random effect, so its drift formula is an intercept.
eamhier_form <- function(arm, bound) {
  mu_rhs <- if (identical(arm, "B")) {
    if (identical(bound, "pg")) {
      frmtmb::bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1 + (1 | s),
                 ndt ~ 1 + (1 | s), bias = 0.5)
    } else {
      frmtmb::bf(rt | dec(upper) ~ 1, bs ~ 1 + (1 | s),
                 ndt ~ 1 + (1 | s), bias = 0.5)
    }
  } else {
    if (identical(bound, "pg")) {
      frmtmb::bf(rt | dec(upper) + ndt_group(s) ~ cond + (1 | s),
                 bs ~ 1 + (1 | s), ndt ~ 1 + (1 | s), bias = 0.5)
    } else {
      frmtmb::bf(rt | dec(upper) ~ cond + (1 | s), bs ~ 1 + (1 | s),
                 ndt ~ 1 + (1 | s), bias = 0.5)
    }
  }
  mu_rhs
}

eamhier_family <- function(arm) {
  if (identical(arm, "C")) {
    frmtmb.eam::wiener(variability = "sv")
  } else {
    frmtmb.eam::wiener()
  }
}

# ------------------------------------------------------- the recording

# The peak process working set in megabytes, read from the operating
# system, because R's gc() cannot see RTMB's tape and the tape is where
# a hierarchical Wiener spends its memory.
eamhier_peak_ws <- function() {
  cmd <- sprintf("(Get-Process -Id %d).PeakWorkingSet64", Sys.getpid())
  out <- tryCatch(suppressWarnings(system2(
    "powershell", c("-NoProfile", "-Command", cmd), stdout = TRUE)),
    error = function(e) NA_character_)
  suppressWarnings(as.numeric(trimws(out[length(out)]))) / 1024^2
}

# One coefficient's Wald interval, by the name confint() gives it. NA
# when the row is absent, so that a missing coefficient reads as missing
# rather than as a silent zero.
eamhier_ci <- function(ci, key) {
  j <- which(rownames(ci) == key)
  if (!length(j)) return(c(NA_real_, NA_real_))
  as.numeric(ci[j[1L], 1:2])
}

eamhier_fmt <- function(v) {
  if (is.numeric(v)) formatC(v, digits = 10, format = "g") else {
    as.character(v)
  }
}

# key=value pairs on one tab-separated line, one file per replicate, so
# that parallel workers never write to the same file.
eamhier_record <- function(path, ...) {
  vals <- list(...)
  line <- paste(paste0(names(vals), "=",
                       vapply(vals, eamhier_fmt, character(1))),
                collapse = "\t")
  cat(line, "\n", sep = "", file = path)
  message("EAMHIER ", line)
  invisible(line)
}
