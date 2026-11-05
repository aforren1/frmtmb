# A mutant to show the pin's cause: read the level's prior from the block density at
# t times a unit vector, t from REV_T, so the difference from the normalizer does
# not round. make: builds the source copy. run: fits on it.
a <- commandArgs(TRUE)
WT <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp"
if (a[[1]] == "make") {
  dst <- "C:/Users/adf44/source/r/rev-impsamp-src/fr1-fix_t"
  unlink(dst, recursive = TRUE)
  dir.create(file.path(dst, "R"), recursive = TRUE)
  dir.create("C:/Users/adf44/source/r/rev-impsamp-lib/fix_t", recursive = TRUE, showWarnings = FALSE)
  for (f in c("DESCRIPTION", "NAMESPACE")) file.copy(file.path(WT, f), dst)
  file.copy(list.files(file.path(WT, "R"), full.names = TRUE), file.path(dst, "R"))
  for (d in c("inst", "data", "man")) file.copy(file.path(WT, d), dst, recursive = TRUE)
  p <- file.path(dst, "R", "importance.R")
  s <- paste0(paste(readLines(p), collapse = "\n"), "\n")
  ed <- list(
    c("  unit <- lapply(seq_len(q), function(j) replace(zero, j, 1))\n",
      "  t0 <- as.numeric(Sys.getenv(\"REV_T\", \"1\"))\n  unit <- lapply(seq_len(q), function(j) replace(zero, j, t0))\n"),
    c("replace(zero, c(j, k), 1)", "replace(zero, c(j, k), t0)"),
    c("dg <- lapply(unit, function(e) bfn(e, theta, bk1) - c0)",
      "dg <- lapply(unit, function(e) (bfn(e, theta, bk1) - c0) / t0^2)"),
    c("value = bfn(pair[[paste0(j, \".\", k)]], theta, bk1) - c0 -",
      "value = (bfn(pair[[paste0(j, \".\", k)]], theta, bk1) - c0) / t0^2 -"))
  for (e in ed) {
    if (!grepl(e[1], s, fixed = TRUE)) stop("no site: ", e[1])
    s <- sub(e[1], e[2], s, fixed = TRUE)
  }
  writeBin(charToRaw(s), p)
  cat("made\n")
} else {
  .libPaths(c("C:/Users/adf44/source/r/rev-impsamp-lib/fix_t", "C:/Users/adf44/source/r/wt-impsamp-lib",
              "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressPackageStartupMessages(library(frmtmb))
  cat(find.package("frmtmb"), "REV_T", Sys.getenv("REV_T"), "\n")
  for (k in c(6, 7, 9)) {
    Sys.setenv(REV_T = format(10^k, scientific = TRUE))
    res <- character(0)
    for (s in 1:8) {
      set.seed(s)
      g <- factor(rep(seq_len(15L), each = 10))
      z <- rnorm(150)
      u <- rnorm(15, 0, 0.7)[g]
      d <- data.frame(y = (0.3 + 0.5 * z + u + rnorm(150, 0, 0.5)) * 10^k, z, g)
      f <- suppressWarnings(tryCatch(frm(bf(y ~ z + (1 | g)), family = gaussian(), data = d, importance = 64L), error = function(e) "STOP"))
      l <- suppressWarnings(frm(bf(y ~ z + (1 | g)), family = gaussian(), data = d))
      res <- c(res, if (is.character(f)) "STOP" else format(as.numeric(logLik(f)) - as.numeric(logLik(l)), digits = 2))
    }
    cat("gaussian y x 1e", k, ":", res, "\n")
    set.seed(1)
    g <- factor(rep(seq_len(15L), each = 10))
    z <- rnorm(150)
    u <- rnorm(15, 0, 0.7)[g]
    d <- data.frame(y = (0.3 + 0.5 * z + u + 0.5 * rt(150, 5)) * 10^k, z, g)
    f <- suppressWarnings(tryCatch(frm(bf(y ~ z + (1 | g), nu = 4), family = student(), data = d, importance = 200L), error = function(e) conditionMessage(e)))
    d$y <- d$y / 10^k
    r <- suppressWarnings(frm(bf(y ~ z + (1 | g), nu = 4), family = student(), data = d, importance = 200L))
    cat("student y x 1e", k, ":", if (is.character(f)) substr(f, 1, 120) else
      sprintf("fit, ll - unit fit %+.4f mcse %.4f", as.numeric(logLik(f)) - (as.numeric(logLik(r)) - 150 * k * log(10)), f$importance$mcse), "\n")
  }
}
