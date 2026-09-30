# Lane sampfix, script 02: where does "no more scalars to read" come from
# on a one-parameter model? Tries the init routes on the reference build.
#
#   Rscript dev/sampfix-02-onepar.R <lane|ref>

arm <- commandArgs(trailingOnly = TRUE)[1L]
LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
cat("ARM ", arm, "\n", sep = "")

set.seed(1212L)
dd <- data.frame(x = rnorm(30))
dd$y <- rpois(30, exp(0.4 * dd$x))
fx <- frm(bf(y ~ 0 + x), family = poisson(), data = dd)
cat("outer parameters: ", length(fx$obj$par), "\n", sep = "")

run <- function(lab, ...) {
  r <- tryCatch(suppressWarnings(suppressMessages(
    frm_sample(fx, chains = 2, iter = 300, refresh = 0, seed = 3, ...))),
    error = function(e) e)
  if (inherits(r, "error")) {
    cat("[", lab, "] ERROR: ", substr(conditionMessage(r), 1, 90), "\n",
        sep = "")
  } else {
    cat("[", lab, "] OK ", paste(dim(r$draws), collapse = "x"), " mean b_x ",
        format(mean(r$draws[, "b_x"]), digits = 6), "\n", sep = "")
  }
}
run("default init (mode)")
run("init = random", init = "random")
run("init_jitter = 0", init_jitter = 0)
# the direct tmbstan call with a length-1 numeric init and with the same
# value carried as a one-dimensional array
obj <- fx$obj
for (lab in c("numeric", "array")) {
  v <- as.numeric(obj$env$last.par.best)
  if (lab == "array") v <- array(v, dim = length(v))
  r <- tryCatch(suppressWarnings(capture.output(
    sf <- tmbstan::tmbstan(obj, chains = 1, iter = 200, refresh = 0,
                           seed = 3, init = list(v)))),
    error = function(e) e)
  ok <- exists("sf") && length(sf@sim) && length(sf@sim$samples)
  cat("[tmbstan init ", lab, "] ", if (ok) "OK" else "NO DRAWS", "\n",
      sep = "")
  if (exists("sf")) rm(sf)
}
cat("DONE\n")
