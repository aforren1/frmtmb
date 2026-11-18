# Lane ordstarts, build round: find a K = 3 ordinal mixture design whose
# default fits converge, for the K >= 3 rule of the swapped start. Data:
# three latent classes (probabilities 0.4, 0.35, 0.25), slopes and shifts
# per class, logistic noise, the cut points of dev/optima-csmix.R; "sep"
# scales the class slopes and shifts. Per seed: the fit with starts = 0,
# its code, start class and elapsed seconds.
#   Rscript dev/ordstarts-k3probe.R <n> <sep> <re: none|int> <seeds> <out>
args <- commandArgs(trailingOnly = TRUE)
nobs <- as.integer(args[1])
sep <- as.numeric(args[2])
re <- args[3]
seeds <- eval(parse(text = args[4]))
out <- args[5]
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
mk3 <- function(seed, G = 50) {
  set.seed(seed)
  x <- rnorm(nobs)
  z <- rnorm(nobs)
  cls <- sample.int(3L, nobs, replace = TRUE, prob = c(0.4, 0.35, 0.25))
  lat <- sep * (c(1.5, -0.8, 0.3)[cls] * x + c(1, -1, 0)[cls]) +
    rlogis(nobs)
  g <- rep_len(seq_len(G), nobs)
  if (re == "int") lat <- lat + rnorm(G, 0, 0.5)[g]
  y <- 1L + (lat > -1 + 0.1 * z) + (lat > 0) + (lat > 1 - 0.1 * z)
  data.frame(x = x, z = z, y = y, g = factor(g))
}
form <- if (re == "int") bf(y ~ x + z + (1 | g)) else bf(y ~ x + z)
q <- function(e) suppressWarnings(suppressMessages(e))
for (s in seeds) {
  d <- mk3(s)
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(q(frm(form, family = mixture(cumulative(), sratio(), acat()),
                      data = d, control = frmtmb_control(starts = 0))),
                error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  row <- if (inherits(f, "error")) {
    data.frame(n = nobs, sep = sep, re = re, seed = s, code = NA, cls = "err",
               ll = NA, el = el)
  } else {
    cl <- ns$mixture_ord_start_class(f, names(f$spec$responses)[1])
    data.frame(n = nobs, sep = sep, re = re, seed = s,
               code = f$opt$convergence,
               cls = if (length(cl$step)) "step" else
                 if (length(cl$boundary)) "flat" else "proper",
               ll = as.numeric(logLik(f)), el = el)
  }
  utils::write.table(row, out, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = !file.exists(out), append = file.exists(out))
  print(row)
}
