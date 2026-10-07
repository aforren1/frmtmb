# Reviewer of the 0.69.0 consolidation: the fuzz tier's generator on
# many seeds of the family of the failing spec (fuzz seed 20379118:
# Gamma, weights, ar1, mo(xo) * z, shape ~ 1 + z, REML), one fit per
# seed, to count fits that stop with an error where 0.68.1 returns one.
#
#   Rscript dev/relrev069-fuzzsweep.R <lib or "base"> <base lib> <variant>
#           <k from> <k to> <out tsv>
# Seeds are 20260901 + k * 977, the plan's own spacing.
# <variant>: "spec" (the failing spec), "mo" (mo(xo) in place of
# mo(xo) * z), "ml" (the spec under ML), "gauss" (gaussian in place of
# Gamma; sigma ~ 1 + z).
# <fix>: optional 7th argument "fix" patches optimize_obj() so a
# restart that raises keeps the run before it (the reviewer's candidate
# fix), to count what it changes.
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(if (!identical(a[1], "base")) a[1], a[2],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
root <- "C:/Users/adf44/source/r/frmtmb-wt-release"
env <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(root, "tests/testthat/helper-fuzz.R"), envir = env)
variant <- a[3]
ks <- seq(as.integer(a[4]), as.integer(a[5]))
out <- a[6]
fix <- length(a) >= 7 && identical(a[7], "fix")
if (fix) {
  ns <- asNamespace("frmtmb")
  oo <- get("optimize_obj", ns)
  src <- deparse(oo)
  i <- grep("opt2 <- run(opt$par)", src, fixed = TRUE)
  stopifnot(length(i) == 1L)
  src[i] <- sub("opt2 <- run(opt$par)",
                "opt2 <- tryCatch(run(opt$par), error = function(e) NULL); if (is.null(opt2)) break",
                src[i], fixed = TRUE)
  oo2 <- eval(parse(text = src))
  environment(oo2) <- ns
  unlockBinding("optimize_obj", ns)
  assign("optimize_obj", oo2, envir = ns)
}
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "| R:", R.home(), "| threads:",
    Sys.getenv("OPENBLAS_NUM_THREADS"), "| variant:", variant,
    "| fix:", fix, "\n")
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
