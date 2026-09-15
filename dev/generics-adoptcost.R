# What does the active binding cost, and where?
#
# The round-1 mechanism paid ONCE PER SESSION at `.onLoad`. The active
# binding pays ONCE PER ACCESS of the name, so the cost model changed
# and has to be re-measured rather than carried over.
#
# CHECK THE INSTRUMENT.  Three things this script does because the
# first spelling of it got them wrong:
#
#  * `proc.time()` ticks at 10 to 20 ms here, so every block is grown
#    past 1.2 s and the statistic is a MINIMUM.
#  * The arms are INTERLEAVED, round by round, rather than run one
#    after another. Run sequentially, the ratio of the binding lookup
#    to one `fixef()` call came out 0.73 in one run and 0.36 in the
#    next on identical code, because `fixef()` drifted from 12 to 22
#    us between arms. A ratio to a drifting arm is not a measurement.
#  * A CONTROL built from the same shape of work, an ordinary binding
#    lookup, is reported beside it, and the headline is the absolute
#    microseconds rather than a ratio to anything that moves.
#
# The round-1 version also measured a defect it had created:
# `setHook()` appends, so timing `frm_adopt_generics()` in a loop timed
# a hook list that grew with the repetition count. There are no hooks.
#
#   Rscript dev/generics-adoptcost.R
LIB <- "C:/Users/adf44/source/r/generics-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))

own <- get("frm_generic_owners", envir = asNamespace("frmtmb"))
pk <- as.environment("package:frmtmb")

block <- function(fun, reps) {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(reps)) fun()
  proc.time()[["elapsed"]] - t0
}
calibrate <- function(fun) {
  reps <- 1024L
  repeat {
    if (block(fun, reps) >= 1.2 || reps > 8e6) break
    reps <- reps * 4L
  }
  reps
}
# one round of every arm, then the next round: machine drift then
# moves all the arms together instead of one of them
run_arms <- function(arms, rounds = 5L) {
  reps <- vapply(arms, calibrate, 1)
  best <- rep(Inf, length(arms))
  for (r in seq_len(rounds)) {
    for (j in seq_along(arms)) {
      best[j] <- min(best[j], block(arms[[j]], reps[j]))
    }
  }
  out <- best / reps * 1e6
  for (j in seq_along(arms)) {
    cat(sprintf("%-42s %9.0f reps %8.2f us\n", names(arms)[j], reps[j],
                out[j]))
  }
  stats::setNames(out, names(arms))
}

tickof <- function(f, want = 200L, budget = 5) {
  d <- numeric(0)
  t0 <- Sys.time()
  repeat {
    v <- replicate(20000, f())
    e <- diff(v)
    d <- c(d, e[e > 0])
    if (length(d) >= want) break
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > budget) {
      break
    }
  }
  if (!length(d)) return(NA_real_)
  min(d)
}

set.seed(1)
dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
dd$y <- rnorm(200, 1 + 0.5 * dd$x, 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

cat("```\n")
cat("== active-binding cost, dev/generics-adoptcost.R ==\n")
cat(sprintf("proc.time() tick here: %.4f s\n",
            tickof(function() proc.time()[["elapsed"]])))
cat(sprintf("names in the table %d, distinct owners %d\n\n",
            length(own), length(unique(unlist(own)))))

cat("-- NO owner loaded: the fallback path, arms interleaved --\n")
a <- run_arms(list(
  `active binding, get('loo')` = function() get("loo", envir = pk),
  `active binding, get('ngrps'), 2 owners` =
    function() get("ngrps", envir = pk),
  `CONTROL get('frm'), ordinary binding` =
    function() get("frm", envir = pk),
  `CONTROL identity(1)` = function() identity(1),
  `fixef(fit), a real accessor call` = function() fixef(fit)))

suppressMessages(loadNamespace("posterior"))
suppressMessages(loadNamespace("loo"))
cat("\n-- posterior and loo loaded: the adopted path --\n")
b <- run_arms(list(
  `active binding, get('loo')` = function() get("loo", envir = pk),
  `active binding, get('as_draws_df')` =
    function() get("as_draws_df", envir = pk),
  `CONTROL get('frm'), ordinary binding` =
    function() get("frm", envir = pk),
  `fixef(fit), a real accessor call` = function() fixef(fit)))

cat(sprintf("\nfallback path %.2f us, adopted path %.2f us\n",
            a[["active binding, get('loo')"]],
            b[["active binding, get('loo')"]]))
cat(sprintf("adopted path is %.1fx an ordinary binding lookup (%.2f us)\n",
            b[["active binding, get('loo')"]] /
              b[["CONTROL get('frm'), ordinary binding"]],
            b[["CONTROL get('frm'), ordinary binding"]]))
cat(sprintf("and %.2f of one fixef() call (%.2f us), same round\n",
            b[["active binding, get('loo')"]] /
              b[["fixef(fit), a real accessor call"]],
            b[["fixef(fit), a real accessor call"]]))
cat("```\n")
