# Reviewer of lane ordstarts, build round: own reproduction of the
# before/after table on two designs (questions 5 and 8).
#   Rscript ordstarts-build-review-rep.R <base|lane> <design> <seeds> <out.tsv>
# design: cs1000 = mixture(cumulative(), sratio()), (1 | g), n 1000, 50 groups
#         cc300  = mixture(cumulative(), cumulative()), (1 | g), n 300, 20 groups
# base (rellib-r22, 0.83.0): arms A (default) and B (starts = 10).
# lane (wt-ordstarts-lib): arms A0 (starts = 0) and D (default).
# My own data generator, written from the findings' description.
args <- commandArgs(trailingOnly = TRUE)
which_lib <- args[1]
design <- args[2]
seeds <- eval(parse(text = args[3]))
out <- args[4]
libs <- c(if (which_lib == "lane") "C:/Users/adf44/source/r/wt-ordstarts-lib",
          "C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
ns <- asNamespace("frmtmb")
n <- if (design == "cs1000") 1000L else 300L
G <- if (design == "cs1000") 50L else 20L
fam <- function() {
  if (design == "cs1000") mixture(cumulative(), sratio()) else
    mixture(cumulative(), cumulative())
}
gen <- function(seed) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  g <- rep_len(seq_len(G), n)
  u0 <- rnorm(G, 0, 0.5)
  lat <- lat + u0[g]
  y <- 1L + (lat > -1 + 0.1 * z) + (lat > 0) + (lat > 1 - 0.1 * z)
  data.frame(y = y, x = x, z = z, g = factor(g))
}
grab <- function(expr) {
  w <- character()
  mc <- character()
  t0 <- proc.time()[["elapsed"]]
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) {
      mc <<- c(mc, class(x)[1])
      invokeRestart("muffleMessage")
    }), error = function(e) e)
  list(v = v, w = w, mc = mc, el = proc.time()[["elapsed"]] - t0)
}
arms <- if (which_lib == "base") {
  list(A = frmtmb_control(), B = frmtmb_control(starts = 10))
} else {
  list(A0 = frmtmb_control(starts = 0), D = frmtmb_control())
}
first <- !file.exists(out)
for (seed in seeds) {
  d <- gen(seed)
  for (a in names(arms)) {
    r <- grab(frm(bf(y ~ x + z + (1 | g)), family = fam(), data = d,
                  control = arms[[a]]))
    f <- r$v
    row <- data.frame(design = design, seed = seed, arm = a, el = r$el,
                      err = inherits(f, "error"), ll = NA, code = NA,
                      nw = NA, snw = NA, degen = NA, lost = NA,
                      maxtse = NA, swapped = NA, swaps = NA, starts = NA,
                      swapmsg = NA)
    if (!inherits(f, "error")) {
      s <- grab(summary(f))
      allw <- c(r$w, s$w)
      fe <- tryCatch(fixef(f), error = function(e) NULL)
      tse <- if (is.matrix(fe)) {
        fe[grepl("Intercept", rownames(fe)), "Est.Error"]
      } else NA
      ms <- f$opt$mixture_starts
      row$ll <- as.numeric(logLik(f))
      row$code <- f$opt$convergence
      row$nw <- length(r$w)
      row$snw <- length(s$w)
      row$degen <- any(grepl("degenerate boundary", allw, fixed = TRUE))
      row$lost <- any(grepl("Standard errors are not available", allw,
                            fixed = TRUE))
      row$maxtse <- suppressWarnings(max(tse, na.rm = TRUE))
      row$swapped <- if (is.null(ms)) 0 else if ("swapped" %in% names(ms)) ms[["swapped"]] else 0
      row$swaps <- if (is.null(ms)) 0 else if ("swaps" %in% names(ms)) ms[["swaps"]] else 0
      row$starts <- if (is.null(ms)) 0 else ms[["starts"]]
      row$swapmsg <- any(r$mc == "frmtmb_mixture_swap")
    }
    utils::write.table(row, out, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = first, append = !first)
    first <- FALSE
  }
  cat(design, which_lib, "seed", seed, "done\n")
}
cat("DONE rep\n")
