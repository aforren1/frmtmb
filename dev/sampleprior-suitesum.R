# Lane sampleprior: per-package totals of the suite runs, from the
# RESULT lines dev/sampleprior-par.sh collected, for the findings.
#   Rscript dev/sampleprior-suitesum.R > dev/sampleprior-log/suites.txt
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-suite"
runs <- if (length(commandArgs(TRUE))) {
  stats::setNames(commandArgs(TRUE)[-1], commandArgs(TRUE)[1])
} else c(`final-1` = "core (ungated) and frmtmb.sample (gated), final build",
          `ext-1` = "the six other extensions (gated), on the lane core",
          `coregated-1` = "core's gated brms files that read priors")
for (r in names(runs)) {
  s <- readLines(file.path(wt, paste0(r, ".summary")))
  pkg <- sub("--.*$", "", s)
  num <- function(k) {
    as.integer(sub(paste0(".* ", k, " ([0-9]+).*"), "\\1", s))
  }
  lane_lib <- grepl("wt-sampleprior-lib/frmtmb ", s, fixed = TRUE)
  d <- data.frame(pkg, files = 1L, pass = num("pass"), fail = num("fail"),
                  error = num("error"), skip = num("skip"),
                  warn = num("warn"), lane_core = lane_lib)
  if (identical(r, "ext-1")) d <- d[d$pkg != "frmtmb.sample", ]
  a <- stats::aggregate(cbind(files, pass, fail, error, skip, warn,
                              lane_core) ~ pkg, data = d, FUN = sum)
  cat(sprintf("%s: %s\n", r, runs[[r]]))
  for (i in seq_len(nrow(a))) {
    cat(sprintf("  %-16s files %3d pass %5d fail %d error %d skip %3d warn %d (lane core in %d of %d)\n",
                a$pkg[i], a$files[i], a$pass[i], a$fail[i], a$error[i],
                a$skip[i], a$warn[i], a$lane_core[i], a$files[i]))
  }
}
