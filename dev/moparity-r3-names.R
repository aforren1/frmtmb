# Punch round 2, r3: data columns and nonlinear parameters NAMED like a
# special (mo, mi, ar, ...) in each formula position the parser reads.
# Prints one line per name and formula: OK with the logLik, or the error.
#   Rscript dev/moparity-r3-names.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args)) args[[1]] else
  "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
nms <- c("mo", "mi", "me", "mmc", "mm", "cs", "gp", "s", "t2", "te", "ar",
         "ma", "arma", "cosy", "unstr", "car", "sar", "fcor", "cens",
         "trunc", "se", "rate", "trials", "dec", "thres", "sm", "rr", "us",
         "diag", "ou", "toep", "cv")
set.seed(2)
n <- 200
base <- data.frame(z = rnorm(n), g = factor(sample(1:15, n, TRUE)))
fails <- 0L
for (nm in nms) {
  d <- base
  d[[nm]] <- rnorm(n)
  d$y <- 1 + 0.3 * d[[nm]] * d$z + rnorm(15)[d$g] + rnorm(n)
  d$yp <- exp(0.5 + 0.4 * d$z) + rnorm(n, 0, 0.2)
  fos <- list(
    str2lang(sprintf("y ~ %s * z", nm)),
    str2lang(sprintf("y ~ %s:z", nm)),
    str2lang(sprintf("y ~ %s + (1 + %s || g)", nm, nm)),
    str2lang(sprintf("y ~ %s + (1 + %s | g)", nm, nm)))
  for (fo in fos) {
    r <- tryCatch(suppressWarnings(suppressMessages(
      frm(bf(eval(fo)), family = gaussian(), data = d))),
      error = function(e) e)
    ok <- !inherits(r, "error")
    if (!ok) fails <- fails + 1L
    cat(sprintf("%-6s %-28s %s\n", nm, deparse1(fo),
                if (ok) paste("OK", format(as.numeric(logLik(r)),
                                           digits = 10))
                else paste("ERROR", substr(conditionMessage(r), 1, 110))))
  }
  # the name as a nonlinear parameter
  nlf <- tryCatch(suppressWarnings(suppressMessages(frm(
    bf(as.formula(sprintf("yp ~ %s + exp(b * z)", nm)),
       as.formula(sprintf("%s ~ 1", nm)), b ~ 1, nl = TRUE),
    family = gaussian(), data = d[setdiff(names(d), nm)],
    prior = set_prior("normal(0, 5)", nlpar = "b")))),
    error = function(e) e)
  ok <- !inherits(nlf, "error")
  if (!ok) fails <- fails + 1L
  cat(sprintf("%-6s %-28s %s\n", nm, "nlpar",
              if (ok) "OK" else paste("ERROR",
                                      substr(conditionMessage(nlf), 1, 110))))
}
cat("failures:", fails, "\n")
