# Lane optrap, consolidation fixes: punch round 2's newton_polish()
# (build 27: no gain gate, never with random effects), read from the
# saved copy of that build's R/fit.R (dev/optrap-c3-base27-fit.R) and
# patched into this process's frmtmb namespace, so the consolidation
# tests can be seen to fail on it with everything else of the final
# build. Sourced after library(frmtmb); dev/optrap-testfile.R does so
# with OPTRAP_ROUND2=1.
local({
  ns <- asNamespace("frmtmb")
  src <- readLines("dev/optrap-c3-base27-fit.R", warn = FALSE)
  s0 <- grep("^newton_polish <- function", src)
  e0 <- s0 - 1L + which(src[s0:length(src)] == "}")[1]
  old <- eval(parse(text = src[s0:e0])[[1]][[3]], envir = ns)
  wrap <- function(obj, opt, bounds, par_units = NULL, control = NULL,
                   ill = FALSE) {
    old(obj, opt, bounds, par_units, control)
  }
  environment(wrap) <- environment()
  utils::assignInNamespace("newton_polish", wrap, "frmtmb")
  cat("newton_polish: punch round 2's (build 27)\n")
})
