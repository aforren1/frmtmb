# Lane optrap, consolidation fixes: the Newton steps' gain gate (keep a
# step only if it gains at least G of its prediction) at another G,
# patched in this process, on the negative binomial check (200 seeds,
# dev/optrap-c3-negbin.R) or the raw polynomials (dev/optrap-c3-repoly.R
# degrees 6 to 8, with or without (1 | g)).
#   OPTRAP_GATE=<G> Rscript dev/optrap-c3-gate.R negbin <lib> <label>
#   OPTRAP_GATE=<G> Rscript dev/optrap-c3-gate.R repoly <lib> <label> <re>
a <- commandArgs(TRUE)
gate <- as.numeric(Sys.getenv("OPTRAP_GATE", "0.25"))
.libPaths(c(a[2], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
b <- deparse(get("newton_polish", ns))
i <- grep("0.25 * gain", b, fixed = TRUE)
stopifnot(length(i) == 1L)
b[i] <- sub("0.25 * gain", paste0(gate, " * gain"), b[i], fixed = TRUE)
# OPTRAP_STRICT=1: undo the steps whenever they end where the Hessian has
# no Cholesky factor, whatever the gradient there
if (nzchar(Sys.getenv("OPTRAP_STRICT"))) {
  j <- grep("g1 > (control$grad_tol", b, fixed = TRUE)
  stopifnot(length(j) == 1L)
  b[j] <- "moved <- FALSE"
  cat("strict undo\n")
}
f <- eval(parse(text = b))
environment(f) <- ns
utils::assignInNamespace("newton_polish", f, "frmtmb")
cat("gain gate patched to", gate, "\n")
args <- if (a[1] == "negbin") c(a[2], a[3]) else c(a[2], a[3], "6:8", a[4])
commandArgs <- function(trailingOnly = FALSE) args
source(if (a[1] == "negbin") "dev/optrap-c3-negbin.R" else
  "dev/optrap-c3-repoly.R")
