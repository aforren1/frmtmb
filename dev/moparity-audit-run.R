# The per-scale floor's audit (dev/moparity-findings.md, item 5): run one
# test file with RTMB::MakeADFun() wrapped so that every objective,
# gradient and Hessian call on an object with a `theta` component
# records the smallest theta it was evaluated at. sd_floored() equals
# exp() bit for bit at a log sd of -308.45 or more, so a file whose
# smallest theta stays above that evaluated every changed block at
# values the floor leaves alone.
#   Rscript dev/moparity-audit-run.R <package> <test file> <out file>
local({
LIB <- "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER"))) Sys.setenv(R_MAKEVARS_USER = mk)
a <- commandArgs(trailingOnly = TRUE)
p <- a[1]
f <- a[2]
out <- a[3]
aud <- new.env()
aud$min <- Inf
aud$n <- 0
orig <- RTMB::MakeADFun
wrapper <- function(...) {
  obj <- orig(...)
  pn <- names(obj$par)
  th <- which(pn == "theta")
  if (!length(th)) return(obj)
  for (nm in c("fn", "gr", "he")) {
    g <- obj[[nm]]
    if (!is.function(g)) next
    obj[[nm]] <- local({
      g0 <- g
      function(...) {
        args <- list(...)
        x <- if (length(args)) args[[1]] else NULL
        if (is.numeric(x) && length(x) == length(pn)) {
          aud$min <- min(aud$min, x[th])
          aud$n <- aud$n + 1
        }
        g0(...)
      }
    })
  }
  obj
}
utils::assignInNamespace("MakeADFun", wrapper, "RTMB")
suppressMessages(library(testthat))
suppressMessages(library(p, character.only = TRUE))
res <- tryCatch(as.data.frame(test_file(f, package = p,
                                        env = testthat::test_env(p),
                                        reporter = "silent")),
                error = function(e) NULL)
cat(sprintf("AUDIT %s %s calls %d min_theta %s fail %s\n", p, basename(f),
            aud$n, format(aud$min, digits = 10),
            if (is.null(res)) "LOADERROR" else sum(res$failed)),
    file = out)
})
