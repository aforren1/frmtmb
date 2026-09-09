# The Phase 0 measurement tier of dev/extension-gaps-plan.md.
#
# One fit per package at the design that plan's "Realistic scale, per
# field" table names, timed and written to a file, so that the cost of a
# realistic model is a number in the repository and not a guess. The
# result of a run is `dev/scale-findings.md`.
#
# To run it:
#
#   Sys.setenv(FRMTMB_SCALE_TESTS = "true", NOT_CRAN = "true",
#              FRMTMB_SCALE_OUT = "<a file to append to>")
#   testthat::test_file("tests/testthat/test-scale.R")
#
# The tier skips without the variable, which is why it can sit in the
# ordinary suite: one of these fits takes minutes, and one of them takes
# longer than that.
#
# This file is BYTE-IDENTICAL in every extension that carries a scale
# row. It is copied rather than made into API for the reason
# frmtmb.sample's helper-sampling.R gives for copying
# `expect_vector_equal()`: a test helper that seven suites share is
# cheaper to duplicate than to promise.

# ------------------------------------------------------------- the gate

skip_unless_scale <- function() {
  testthat::skip_on_cran()
  if (!identical(Sys.getenv("FRMTMB_SCALE_TESTS"), "true")) {
    testthat::skip("set FRMTMB_SCALE_TESTS=true to run the scale tier")
  }
}

# A file may hold more than one design. FRMTMB_SCALE_ROW runs one of
# them alone, because a timed fit belongs in a process that has not
# already built a tape of its own.
scale_row_on <- function(row) {
  want <- Sys.getenv("FRMTMB_SCALE_ROW", "")
  if (nzchar(want) && !identical(want, row)) {
    testthat::skip(paste0("FRMTMB_SCALE_ROW is set to '", want, "'"))
  }
  invisible(TRUE)
}

# A STRUCTURAL run: every design shrunk to a size that checks the code
# path without paying for the measurement, so that the tier itself can
# be exercised in a minute. Every row recorded under it carries
# `small=TRUE`, because nothing measured at a toy size is a measurement.
scale_small <- function() {
  identical(Sys.getenv("FRMTMB_SCALE_SMALL"), "true")
}

# ------------------------------------------------------- the instrument

# proc.time() ticks at 10.0 ms on the machine this tier was written on,
# so a gradient measured with it is a handful of ticks and a ratio of
# two such readings is noise. Sys.time() resolves to about 2 us here,
# which is why nothing below reads proc.time().
scale_elapsed <- function(expr) {
  gc(FALSE)
  t0 <- Sys.time()
  force(expr)
  as.numeric(difftime(Sys.time(), t0, units = "secs"))
}

# One gradient evaluation. The batch DOUBLES until the block passes
# `min_block` seconds, so the answer never rests on a few clock ticks,
# and the best of `blocks` blocks is reported, because a block can only
# be made slower by something outside the measurement and never faster.
# `max_calls` is the escape for a gradient that is already seconds long,
# where one call is its own block.
scale_grad <- function(obj, par, min_block = 1.2, blocks = 3L,
                       max_calls = 256L) {
  obj$gr(par)
  n <- 1L
  t <- scale_elapsed(for (i in seq_len(n)) obj$gr(par))
  while (t < min_block && n < max_calls) {
    n <- min(max_calls, 2L * n)
    t <- scale_elapsed(for (i in seq_len(n)) obj$gr(par))
  }
  best <- t / n
  for (b in seq_len(max(0L, blocks - 1L))) {
    best <- min(best, scale_elapsed(for (i in seq_len(n)) obj$gr(par)) / n)
  }
  list(seconds = best, calls = n)
}

# The tape build: the objective time less the frame time, because frame
# assembly is not what a tape build costs.
#
# TWO ROUNDS, interleaved, minimum of each arm. The FIRST call through
# frm() in a process pays one-time costs that belong to neither arm
# (namespace loading, mgcv, method dispatch, the byte-code compiler),
# and a single round of the subtraction was measured NEGATIVE because
# of them on the small coupling designs. The objective object of the
# last round comes back, so the gradient is timed on the tape that was
# just built rather than on another one.
scale_build <- function(form, ..., rounds = 2L) {
  fr <- Inf
  ob <- Inf
  dry <- NULL
  for (r in seq_len(rounds)) {
    fr <- min(fr, scale_elapsed(
      frmtmb::frm(form, ..., dry_run = "frame")))
    ob <- min(ob, scale_elapsed(
      dry <- frmtmb::frm(form, ..., dry_run = "objective")))
  }
  list(frame_s = fr, objective_s = ob, build_s = max(ob - fr, 0),
       dry = dry)
}

# Post-fit calls, timed as a SET rather than one at a time.
#
# `fs` is a named list of zero-argument functions. One call of each is
# made per round, in order, for `rounds` rounds. Interleaved because
# that is what the standing rules require of anything a ratio will be
# taken between, and replicated because a ratio of two single timings
# is not a measurement: this tier shipped one such ratio (12.75x for a
# draw ratio of 10) that did not survive replication, and three
# interleaved rounds gave 10.5.
#
# THREE numbers come back per arm, not one, because for a memoized call
# they are three different questions. `first` is the cold call, which
# is what a user who calls it once pays. `seconds` is the minimum,
# which is the warm call. `spread` is max over min, which is how a
# memoized call announces itself: frm_curve(simultaneous = TRUE) came
# back with a spread of 8.7 and it was the joint-precision solve being
# cached on the fit, not noise.
scale_interleave <- function(fs, rounds = 3L) {
  ts <- matrix(NA_real_, rounds, length(fs),
               dimnames = list(NULL, names(fs)))
  for (r in seq_len(rounds)) {
    for (j in seq_along(fs)) ts[r, j] <- scale_elapsed(fs[[j]]())
  }
  list(seconds = apply(ts, 2L, min),
       first = ts[1L, ],
       spread = apply(ts, 2L, max) / apply(ts, 2L, min),
       rounds = rounds)
}

# The instrument reading itself. Two blocks of the SAME work through the
# same path must report a ratio near one. A control far from one says
# the clock or the machine, and not the model, is what the row measured.
scale_control <- function(obj, par, calls) {
  a <- scale_elapsed(for (i in seq_len(calls)) obj$gr(par))
  b <- scale_elapsed(for (i in seq_len(calls)) obj$gr(par))
  max(a, b) / min(a, b)
}

# Peak R heap since the last reset, in megabytes. It does not see RTMB's
# tape, which lives in C++ memory, so it is a FLOOR on what the fit
# costs and the findings report it as one. The process peak is taken
# from outside, by the runner.
scale_mem_reset <- function() invisible(gc(reset = TRUE))
scale_mem_peak_mb <- function() sum(gc()[, 6L])

# ---------------------------------------------------------- the record

# Recovery is asserted as a z score against the fit's OWN standard
# error, never as an absolute difference. Three absolute tolerances have
# broken on this project's CI already.
scale_z <- function(est, se, truth) abs(est - truth) / se

# What diagnose() found, as one short string, so that a row that did not
# converge still says why.
#
# `unbounded` is core 0.55.0's check for a dpar standing at the edge of
# its link. It is here because the plan's item 2.1 predicts exactly that
# shape of failure for a random effect on a bounded non-decision time,
# and a row that reproduces it should say so rather than leave a reader
# to infer it from a large gradient.
scale_diag <- function(fit) {
  d <- frmtmb::diagnose(fit, quiet = TRUE)
  # a data.frame of offending coefficients, or NULL. Note that it looks
  # at a dpar's FIXED-effect columns, so a variance component standing
  # at the edge of the same link is not what it reports.
  ub <- d$unbounded_dpar
  ub <- if (is.null(ub)) "" else paste(ub$parameter, collapse = "/")
  paste0("conv=", d$convergence, ",maxgrad=",
         formatC(d$max_grad, digits = 3, format = "g"),
         ",pdHess=", isTRUE(d$pdHess),
         ",nbadse=", length(d$bad_se),
         ",nflat=", length(d$flat),
         ",unbounded=", if (nzchar(ub)) ub else "none")
}

# One line per row, key=value pairs separated by tabs, appended, so that
# several rows run in several processes build one table.
scale_record <- function(row, ...) {
  vals <- list(...)
  fmt <- function(v) {
    if (is.numeric(v)) formatC(v, digits = 6, format = "g") else {
      as.character(v)
    }
  }
  line <- paste(c(paste0("row=", row),
                  paste0("small=", scale_small()),
                  paste0(names(vals), "=",
                         vapply(vals, fmt, character(1)))),
                collapse = "\t")
  path <- Sys.getenv("FRMTMB_SCALE_OUT", "")
  if (!nzchar(path)) path <- file.path(tempdir(), "frmtmb-scale.tsv")
  cat(line, "\n", sep = "", file = path, append = TRUE)
  message("SCALE ", line)
  invisible(line)
}
