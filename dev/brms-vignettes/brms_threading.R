# brms's "Running brms models with within-chain parallelization"
# (brms 2.23.0, doc/brms_threading.Rmd) against the frmtmb surface.
#
# This vignette is NOT translated. It is about one thing only: splitting
# the evaluation of a Stan log-likelihood across threads inside a single
# NUTS chain. Every line of it is either a Stan compilation option, a
# cmdstanr backend argument, or a benchmark harness for the two. None of
# that has a frmtmb counterpart, and the reason is structural rather
# than a gap waiting to be filled, so the script records the surface and
# the reason instead of attempting a port.
#
# Labels carry the "ML: " prefix so the summarizer's path split stays
# well defined. There is no sampling section: the vignette's subject is
# the sampler's internals, not a model.

source(file.path(Sys.getenv(
  "BV_DIR",
  unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig/dev/brms-vignettes"
), "_harness.R"))
bv_load()
bv_init("brms_threading")

options(mc.cores = 1)

bv("data", "the vignette's subject", {
  cat("brms_threading is a Stan within-chain performance study.\n")
  cat("It runs one model repeatedly under varying thread counts and\n")
  cat("grainsizes and plots the speedup. Nothing in it is a modeling\n")
  cat("question, so nothing in it is translated.\n")
  "not translated"
}, NA_character_, "")

bv("post", "ML: threading(threads, grainsize, static)", {
  threading(4)
}, "MISSING",
"threading() has no frmtmb analog: it configures Stan's reduce_sum and there is no Stan program. Nothing in frmtmb names a replacement because there is no equivalent knob")

bv("post", "ML: brm(..., threads = threading(2))", {
  stop("`threads` is not an argument of frm() or frm_sample()")
}, "MISSING",
"the `threads` argument does not exist on either frmtmb fitting entry point; a user copying the line meets an unused-argument error from frm()")

bv("post", "ML: reduce_sum in the generated Stan code", {
  stop("there is no generated Stan code to add reduce_sum to")
}, "MISSING",
"reduce_sum partitions a Stan log-likelihood sum. frmtmb's likelihood is an R closure taped by RTMB, and RTMB does not expose TMB's OpenMP parallel_accumulator, so the partition point does not exist")

bv("post", "ML: brm(..., backend = 'cmdstanr')", {
  stop("`backend` is not an argument of frm() or frm_sample()")
}, "MISSING",
"there is one backend. frm_sample() goes through tmbstan and rstan, and threading is a cmdstanr-only feature in brms, so the argument has nothing to select")

bv("post", "ML: grainsize tuning", {
  stop("grainsize is a reduce_sum argument")
}, "MISSING", "no reduce_sum means no grainsize; nothing to tune")

# What frmtmb does offer, read off frmtmb_control() and dev/benchmarks.md.
bv("post", "ML: what frmtmb offers instead of threading()", {
  ctl <- names(frmtmb_control())
  stopifnot(!any(grepl("thread|core|parallel|omp", ctl, ignore.case = TRUE)))
  cat("frmtmb_control() arguments:", paste(ctl, collapse = ", "), "\n")
  cat("\nNo parallelism argument of any kind. The measured reasons, from\n")
  cat("dev/benchmarks.md:\n")
  cat(" - optimParallel was benchmarked and rejected: it is slower than\n")
  cat("   the sequential optimizer whenever the gradient is exact, which\n")
  cat("   it always is here.\n")
  cat(" - RTMBp with autopar taping was benchmarked and rejected for\n")
  cat("   now: no gain at all on Cholesky-dominated mixed models, and\n")
  cat("   1.7x to 2x end to end only on large accumulation-dominated\n")
  cat("   GLMMs, from an off-CRAN package coupled to TMB internals.\n")
  cat(" - The speedups that do exist are different in kind:\n")
  cat("   frmtmb_control(profile = TRUE) removes the fixed effects from\n")
  cat("   the outer problem, and better start values cut the iteration\n")
  cat("   count. Both attack the same seconds threading would.\n")
  "recorded"
}, "MISSING",
"frmtmb has no parallelism option: frmtmb_control() carries optimizer, optCtrl, restarts, grad_tol, profile, sparse_x, autoscale, two checks and verbose, and none of them is a thread count")

# The one place concurrency does appear, and it is between fits rather
# than inside one.
bv("post", "ML: frm_sample(cores =) is the only parallel argument", {
  a <- names(formals(frm_sample))
  stopifnot(!"cores" %in% a)
  cat("frm_sample() arguments:", paste(a, collapse = ", "), "\n")
  cat("`cores` reaches rstan through ..., so it parallelizes CHAINS,\n")
  cat("not the inside of one chain, and R/interop.R refuses more than\n")
  cat("one core on Windows because rstan uses PSOCK workers there.\n")
  "recorded"
}, "BEHAVIOR",
"the nearest thing to threading() is between-chain parallelism through frm_sample(cores =), which is the axis brms's vignette explicitly says it is NOT about, and it is disabled on Windows")

bv("post", "ML: fit-level parallelism is the documented substitute", {
  cat("dev/benchmarks.md's conclusion: for repeated fits over one\n")
  cat("dataset (frm_bootstrap, profile, frm_allfit) a persistent\n")
  cat("cluster running WHOLE fits is where the concurrency is. None of\n")
  cat("those functions takes a cores argument today, so the user builds\n")
  cat("the cluster themselves.\n")
  stopifnot(!"cores" %in% names(formals(frm_bootstrap)),
            !"cores" %in% names(formals(frm_allfit)))
  "recorded"
}, "MISSING",
"frm_bootstrap() and frm_allfit() are the embarrassingly parallel entry points and neither takes cores; the parallelism benchmarks.md recommends is left to the caller")

bv_done()
