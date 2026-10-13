# Lane optrap: dev/relrev069-fuzzsweep.R run from this worktree (its
# helper-fuzz.R), without the candidate patch: the lane's build carries
# its own fix. Everything else is that script's.
#   Rscript dev/optrap-fuzzsweep.R <lib or "base"> <base lib> <variant>
#           <k from> <k to> <out tsv>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (!identical(a[1], "base")) a[1], a[2],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
root <- "C:/Users/adf44/source/r/frmtmb-wt-optrap"
env <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(root, "tests/testthat/helper-fuzz.R"), envir = env)
variant <- a[3]
ks <- seq(as.integer(a[4]), as.integer(a[5]))
out <- a[6]
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "| R:", R.home(), "| threads:",
    Sys.getenv("OPENBLAS_NUM_THREADS"), "| variant:", variant, "\n")
sp0 <- list(family = "Gamma", aterm = "weights", re = "ar1",
            special = "mo_int", dpar = "dpar_x", mode = "reml",
            op = "confint")
if (variant == "mo") sp0$special <- "mo"
if (variant == "ml") sp0$mode <- "ml"
if (variant == "gauss") sp0$family <- "gaussian"
cat("seed\tstatus\tcode\tobjective\tnwarn\tsecs\tdetail\n", file = out)
for (k in ks) {
  sp <- sp0
  sp$seed <- 20260901L + k * 977L
  d <- env$fuzz_data(sp)
  if (is.null(d)) {
    cat(sp$seed, "\tnodata\tNA\tNA\tNA\tNA\t\n", sep = "", file = out,
        append = TRUE)
    next
  }
  t0 <- proc.time()[["elapsed"]]
  res <- env$fuzz_fit_one(sp, d)
  el <- proc.time()[["elapsed"]] - t0
  if (is.null(res$value)) {
    msg <- conditionMessage(res$error %||% simpleError("?"))
    cat(sp$seed, "\terror\tNA\tNA\t", length(res$warnings), "\t",
        sprintf("%.1f", el), "\t",
        gsub("[\t\r\n]+", " ", substr(msg, 1, 120)), "\n", sep = "",
        file = out, append = TRUE)
  } else {
    f <- res$value
    cat(sp$seed, "\tok\t", f$opt$convergence, "\t",
        sprintf("%.10g", f$opt$objective), "\t", length(res$warnings),
        "\t", sprintf("%.1f", el), "\t\n", sep = "", file = out,
        append = TRUE)
  }
}
cat("DONE\n")
