# Reviewer's mutation runner for lane bin2port: puts a deliberately
# wrong frmtmb side into a scratch copy of dev/brmsport-bin2-adapters.R
# and runs one generated bin-2 file in tier mode (recorded verdicts and
# value-check outcomes asserted), printing every failing expectation.
# Run with the working directory at a COPY of the lane's tree.
#
#   Rscript <this> <mutation> <topic> [<generated file override>]
#
# Mutations (each a fixed-string substitution that must match once):
#   none      the lane's adapter unchanged
#   catswap   brms's mu dpars mapped onto frmtmb's in reverse order
#   catlpdf   frmtmb's categorical-like row density shifted by 1e-6
#   frechet   frmtmb's frechet mean off by a factor 1.001
#   gev       frmtmb's GEV sigma off by a factor 1.001
#   lncor     the logistic-normal correlations negated
#   lncorT    the logistic-normal correlation matrix filled by columns
lib <- "C:/Users/adf44/source/r/rellib-r9"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
a <- commandArgs(trailingOnly = TRUE)
mut <- a[1]
topic <- a[2]
gen <- if (length(a) >= 3L) a[3] else
  sprintf("dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-%s.R", topic)
subs <- list(
  none = NULL,
  catswap = c("    stats::setNames(b, f_names)\n",
              "    stats::setNames(rev(b), f_names)\n"),
  catlpdf = c("    fam <- fam[[\"defer\"]](Y ~ 1, d)\n  }\n  fam\n}",
              paste0("    fam <- fam[[\"defer\"]](Y ~ 1, d)\n  }\n",
                     "  f0 <- fam$lpdf\n",
                     "  fam$lpdf <- function(...) f0(...) + 1e-6\n",
                     "  fam\n}")),
  frechet = c("list(mu = scale * gamma(1 - 1 / shape), nu = shape)",
              "list(mu = 1.001 * scale * gamma(1 - 1 / shape), nu = shape)"),
  gev = c("brms_g1_d(\"gen_extreme_value\", x, list(mu = mu, sigma = sigma,",
          "brms_g1_d(\"gen_extreme_value\", x, list(mu = mu, sigma = 1.001 * sigma,"),
  lncor = c("C[a, b] <- C[b, a] <- cors[s, k]",
            "C[a, b] <- C[b, a] <- -cors[s, k]"),
  lncorT = c("for (a in 2:D) for (b in seq_len(a - 1L)) {",
             "for (b in 1:(D - 1L)) for (a in (b + 1L):D) {"))
stopifnot(mut %in% c(names(subs), "cumfix"))
src <- paste(readLines("dev/brmsport-bin2-adapters.R"), collapse = "\n")
if (identical(mut, "cumfix")) {
  # fix round 1: the cumulative fix applied in the test session only, so
  # the defect the value check log_lik:462 records is gone
  src <- paste(src, paste(
    "local({",
    "  ns <- asNamespace('frmtmb')",
    "  f <- get('ord_cumulative_logpmf', ns)",
    "  s <- deparse(f)",
    "  p <- 'a <- q(disc * (tau[ym] - eta))'",
    "  stopifnot(sum(grepl(p, s, fixed = TRUE)) == 1L)",
    "  g <- eval(parse(text = sub(p, paste(p, '+ (i1 + iK)'), s,",
    "                             fixed = TRUE)))",
    "  environment(g) <- ns",
    "  unlockBinding('ord_cumulative_logpmf', ns)",
    "  assign('ord_cumulative_logpmf', g, envir = ns)",
    "  lockBinding('ord_cumulative_logpmf', ns)",
    "})", sep = "\n"), sep = "\n")
}
s <- subs[[mut]]
if (!is.null(s)) {
  hits <- gregexpr(s[1], src, fixed = TRUE)[[1]]
  stopifnot(length(hits) == 1L, hits[1] > 0)
  src <- sub(s[1], s[2], src, fixed = TRUE)
}
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
Sys.setenv(R_MAKEVARS_USER = mk, NOT_CRAN = "true",
           FRMTMB_BRMS_FIT_TESTS = "true", FRMTMB_BRMSPORT_PKG = "frmtmb",
           FRMTMB_PORT_ROOT = getwd())
stage <- file.path("dev/bin2port-rev-stage", paste0(mut, "-", topic))
unlink(stage, recursive = TRUE)
dir.create(stage, recursive = TRUE)
file.copy(list.files("tests/testthat", "^helper.*[.]R$", full.names = TRUE),
          stage)
writeLines(src, file.path(stage, "helper-brms-bin2-adapters.R"))
file.copy(gen, file.path(stage, basename(gen)))
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
cat("mutation:", mut, " file:", gen, " lib:",
    dirname(system.file(package = "frmtmb")),
    format(utils::packageVersion("frmtmb")), "\n")
res <- test_file(normalizePath(file.path(stage, basename(gen)),
                               winslash = "/"),
                 package = "frmtmb", env = testthat::test_env("frmtmb"),
                 reporter = "silent", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s %s pass=%d fail=%d err=%d\n", mut, topic,
            sum(df$passed), sum(df$failed), sum(df$error)))
for (t in res) for (e in t$results) {
  if (inherits(e, c("expectation_failure", "expectation_error"))) {
    cat("  FAIL [", t$test, "] ",
        substr(gsub("\\s+", " ", conditionMessage(e)), 1, 260), "\n",
        sep = "")
  }
}
unlink(stage, recursive = TRUE)
