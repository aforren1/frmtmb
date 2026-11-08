# Before and after timings for the rewritten sites that run more than
# once per fit, measured in ONE process so that the two arms see the same
# core type and load (dev/lane-rules.md, "Check the instrument").
#
# The "before" arm is the base source of the same function (frmtmb
# 3cd2a0c8, read from the main checkout's R/ files) evaluated in the lane
# namespace, so both arms call the same helpers. The CONTROL arm is the
# lane's source of the function built the same way. The before/after
# comparison is base against control, two arms built identically; the
# lane arm (the installed function) is shown beside them, and its ratio
# to control measures what the construction itself costs.
#
# Each block repeats a call until it runs at least 1.2 s; arms are
# interleaved and the minimum per-call time over ROUNDS blocks is kept.
#
# Usage: Rscript dev/superassign-bench.R <base source root>

a <- commandArgs(trailingOnly = TRUE)
base_root <- a[1L]
.libPaths(c("C:/Users/adf44/source/r/wt-superassign-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(20261106)
ROUNDS <- 7L

old_def <- function(file, name, root = base_root) {
  ex <- parse(file.path(root, "R", file), keep.source = FALSE)
  for (e in ex) {
    if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
          identical(e[[2L]], as.name(name))) {
      f <- eval(e[[3L]], ns)
      environment(f) <- ns
      # byte-compiled as the installed package is, so that the base and
      # control arms are built the same way
      return(compiler::cmpfun(f))
    }
  }
  stop("no definition of ", name, " in ", file)
}
# the control is built the way the base arm is, from the lane's source
lane_root <- "C:/Users/adf44/source/r/frmtmb-wt-superassign"
lane_copy <- function(name) {
  file <- c(nlminb_trial_fn = "fit.R", frm_bind_generic = "generic-owners.R",
            compat_hand_rules_tbl = "compat.R",
            prior_table = "priors.R")[[name]]
  old_def(file, name, root = lane_root)
}

# per-call seconds of `call(arm)`, minimum over ROUNDS interleaved blocks
bench <- function(label, arms, call) {
  reps <- 1L
  repeat {
    t <- system.time(for (i in seq_len(reps)) call(arms[[1L]]))[["elapsed"]]
    if (t >= 1.2) break
    reps <- reps * 2L
  }
  best <- setNames(rep(Inf, length(arms)), names(arms))
  for (r in seq_len(ROUNDS)) {
    for (k in sample(names(arms))) {
      t <- system.time(for (i in seq_len(reps)) call(arms[[k]]))
      best[[k]] <- min(best[[k]], t[["elapsed"]] / reps)
    }
  }
  cat(sprintf("%s reps=%d\n", label, reps))
  for (k in names(arms)) {
    cat(sprintf("  %-8s %.3e s/call  ratio to lane %.3f\n", k, best[[k]],
                best[[k]] / best[["lane"]]))
  }
  invisible(best)
}

# 1. nlminb_trial_fn(): runs once per objective evaluation the optimizer
#    asks for. The objective here is trivial, so the wrapper's own cost
#    is what is timed.
mk <- list(base = old_def("fit.R", "nlminb_trial_fn"),
           lane = get("nlminb_trial_fn", envir = ns),
           control = lane_copy("nlminb_trial_fn"))
obj <- function(p) sum(p * p)
w <- lapply(mk, function(m) m(obj)$fn)
p0 <- c(0.1, -0.2, 0.3)
bench("nlminb_trial_fn()$fn per call", w, function(f) f(p0))

# 2. frm_bind_generic(): the active binding runs on every access to a
#    borrowed generic. Owner "stats" is loaded, so each access takes
#    the memo branch.
bind <- list(base = old_def("generic-owners.R", "frm_bind_generic"),
             lane = get("frm_bind_generic", envir = ns),
             control = lane_copy("frm_bind_generic"))
envs <- lapply(bind, function(b) {
  e <- new.env()
  b(e, "fitted", "stats", function(...) NULL)
  e
})
stopifnot(identical(envs$base$fitted, stats::fitted),
          identical(envs$lane$fitted, stats::fitted))
bench("active binding access", envs, function(e) e$fitted)

# 3. compat_hand_rules_tbl(): about 600 rows appended one at a time;
#    once per session (cached), timed with the cache cleared.
hand <- list(base = old_def("compat.R", "compat_hand_rules_tbl"),
             lane = get("compat_hand_rules_tbl", envir = ns),
             control = lane_copy("compat_hand_rules_tbl"))
cache <- get("compat_hand_cache", envir = ns)
saved <- cache$tbl
tabs <- lapply(hand, function(h) {
  cache$tbl <- NULL
  h()
})
cat("compat rows:", nrow(tabs$lane), " identical base/lane:",
    identical(tabs$base, tabs$lane), "\n")
bench("compat_hand_rules_tbl() uncached", hand, function(h) {
  cache$tbl <- NULL
  h()
})
cache$tbl <- saved

# 4. prior_table(): the default prior table, built row by row.
dd <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60),
                 g = factor(rep(1:6, 10)), h = factor(rep(1:5, 12)))
fit <- suppressMessages(suppressWarnings(
  frm(bf(y ~ x + z + (1 + x | g) + (1 | h), sigma ~ x),
      family = gaussian(), data = dd)))
pt <- list(base = old_def("priors.R", "prior_table"),
           lane = get("prior_table", envir = ns),
           control = lane_copy("prior_table"))
rows <- lapply(pt, function(f) f(fit$spec, fit$frame, "fit"))
cat("prior rows:", nrow(rows$lane), " identical base/lane:",
    identical(rows$base, rows$lane), "\n")
bench("prior_table()", pt, function(f) f(fit$spec, fit$frame, "fit"))
