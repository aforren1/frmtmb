# Lane optima, item 3: cs() on a cumulative component of an ordinal
# mixture. The release review's construction (dev/relrev-csfail.R): 20
# seeds, n = 300, two-class latent data with no cs() effect, and its
# controls. Each fit is classified as nan_grad (the optimizer's
# "NA/NaN gradient evaluation"), error, nonconv (code not 0), warn (a
# warning other than brms's experimental cs() one) or ok, and its
# logLik is kept.
#   Rscript dev/optima-csmix.R base|lane out.tsv [designs] [seeds]
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
out <- args[2]
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
mk <- function(seed, n = 300, mix = FALSE) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  lat <- if (mix) {
    cls <- rbinom(n, 1, 0.4)
    ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  } else 0.8 * d$x + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
one <- function(expr) {
  w <- character()
  t0 <- proc.time()[[3]]
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }), error = function(e) structure(conditionMessage(e), class = "err"))
  dt <- proc.time()[[3]] - t0
  if (inherits(r, "err")) {
    cls <- if (grepl("NA/NaN gradient", r)) "nan_grad" else "error"
    return(data.frame(class = cls, ll = NA, code = NA, evals = NA,
                      sec = dt, msg = substr(gsub("[\t\n]", " ", r), 1, 200)))
  }
  ow <- w[!grepl("Category specific effects for this family", w,
                 fixed = TRUE)]
  cls <- if (r$opt$convergence != 0) {
    "nonconv"
  } else if (length(ow)) "warn" else "ok"
  data.frame(class = cls, ll = as.numeric(logLik(r)),
             code = r$opt$convergence, evals = r$opt$evals %||% NA,
             sec = dt,
             msg = substr(gsub("[\t\n]", " ", paste(ow, collapse = " | ")),
                          1, 200))
}
designs <- list(
  cum = function(d) frm(y ~ x + cs(z), family = cumulative(), data = d),
  sratio = function(d) frm(y ~ x + cs(z), family = sratio(), data = d),
  mix_cum_sratio = function(d) {
    frm(bf(y ~ x + cs(z)), family = mixture(cumulative(), sratio()),
        data = d)
  },
  mix_cum_cum = function(d) {
    frm(bf(y ~ x + cs(z)), family = mixture(cumulative(), cumulative()),
        data = d)
  },
  mix_sratio_sratio = function(d) {
    frm(bf(y ~ x + cs(z)), family = mixture(sratio(), sratio()), data = d)
  },
  mix_cum_sratio_nocs = function(d) {
    frm(bf(y ~ x + z), family = mixture(cumulative(), sratio()), data = d)
  }
)
sel <- if (length(args) > 2 && nzchar(args[3])) {
  strsplit(args[3], ",")[[1]]
} else names(designs)
seeds <- if (length(args) > 3) eval(parse(text = args[4])) else 1:20
rows <- list()
for (nm in sel) {
  mix <- startsWith(nm, "mix")
  for (s in seeds) {
    r <- one(designs[[nm]](mk(s, mix = mix)))
    rows[[length(rows) + 1L]] <- cbind(design = nm, seed = s, r)
  }
  X <- do.call(rbind, rows)
  tb <- table(X$class[X$design == nm])
  cat(sprintf("%-22s %s\n", nm, paste(names(tb), tb, sep = "=",
                                      collapse = " ")))
}
utils::write.table(do.call(rbind, rows), out, sep = "\t", quote = FALSE,
                   row.names = FALSE)
