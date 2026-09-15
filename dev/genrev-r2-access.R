# genrev round 2: the per-ACCESS cost of the active binding, with the
# instrument checked.  One process per owner state.
#
#  * blocks grown until each takes >= 1.5 s on the elapsed clock, whose
#    tick is measured and printed;
#  * arms interleaved within each of 7 rounds, minimum per arm;
#  * a CONTROL built from the same code that must report the same as
#    frmtmb's own binding: frm_bind_generic() installed, from frmtmb's
#    namespace, into a scratch environment for the same name;
#  * a CONTROL that isolates active-binding overhead: an active binding
#    returning a constant.
a <- commandArgs(trailingOnly = TRUE)
LIB <- a[1]; state <- a[2]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
if (state == "owners") {
  suppressMessages({ loadNamespace("loo"); loadNamespace("posterior") })
}
ns <- asNamespace("frmtmb")
tk <- diff(unique(replicate(200000, proc.time()[["elapsed"]])))
cat(sprintf("state %s  elapsed tick min %.4f s median %.4f s\n",
            state, min(tk), median(tk)))
scratch <- new.env()
get("frm_bind_generic", envir = ns)(scratch, "loo", "loo",
                                     get("loo", envir = ns))
const <- new.env()
fn <- function(x) x
makeActiveBinding("k", function() fn, const)
plain <- new.env(); assign("p", fn, envir = plain)
arms <- list(
  `frmtmb loo (active)` = function() get("loo", envir = ns),
  `frmtmb ngrps (2 owners)` = function() get("ngrps", envir = ns),
  `CONTROL same code, scratch env` = function() get("loo", envir = scratch),
  `CONTROL active constant` = function() get("k", envir = const),
  `CONTROL plain binding` = function() get("p", envir = plain),
  `CONTROL identity` = function() identity(1))
size <- function(f) {
  reps <- 1024L
  repeat {
    el <- system.time(for (i in seq_len(reps)) f())[["elapsed"]]
    if (el >= 1.5) return(reps)
    reps <- reps * 2L
  }
}
reps <- vapply(arms, size, 1L)
best <- setNames(rep(Inf, length(arms)), names(arms))
for (r in 1:7) for (j in sample(seq_along(arms))) {
  el <- system.time(for (i in seq_len(reps[j])) arms[[j]]())[["elapsed"]]
  best[j] <- min(best[j], el / reps[j])
}
for (j in seq_along(arms)) {
  cat(sprintf("  %-32s %9d reps  %6.2f us  (%.0f ticks per block)\n",
              names(arms)[j], reps[j], 1e6 * best[j],
              best[j] * reps[j] / min(tk)))
}
cat(sprintf("frmtmb loo / plain binding = %.2f ; same-code control / frmtmb loo = %.3f\n",
            best[1] / best[5], best[3] / best[1]))
cat("GENREVDONE\n")
