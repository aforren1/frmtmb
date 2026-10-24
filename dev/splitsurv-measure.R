# Guardrail 3 of dev/core-boundary.md, measured BEFORE the moves, for
# (a) curve inference, (b) Royston-Parmar and (c) cox():
# lines of code in R/, test time and its share of the source suite.
# Commits, dependencies and ported brms rows are in the findings, from
# the commands named there.
#
#   Rscript dev/splitsurv-measure.R
#
# Sources: the base commit 19b162b9 (git show), and the per-test times
# of the before run on rellib-r11 (dev/splitsurv-before-out/*.tsv,
# written by dev/splitsurv-par.sh with 8 processes at once on a shared
# machine, so a time is an order of magnitude, not a benchmark).
base <- "19b162b9"
show <- function(path) {
  system2("git", c("show", paste0(base, ":", path)), stdout = TRUE)
}
loc <- function(x) {
  code <- x[!grepl("^\\s*(#|$)", x)]
  c(lines = length(x), code = length(code))
}
sp <- "extensions/frmtmb.spline/R/"
cat("== lines in R/ (lines, of which code: not blank, not a comment)\n")
curve <- unlist(lapply(paste0(sp, c("curve.R", "curve-cov.R",
                                    "curve-deriv.R", "curve-feature.R")),
                       show))
rp <- unlist(lapply(paste0(sp, c("royston-parmar.R", "rp-basis.R",
                                 "rp-check.R")), show))
zzz <- show(paste0(sp, "zzz.R"))
i_cur <- grep("## ---- curve inference", zzz)
i_rpf <- grep("r[(]\"rp_floored\", \"cens[(][)]\"", zzz)
i_end <- grep("^  b[$]rules[(][)]", zzz)
i_fam <- grep("## ---- the family", zzz)
rp_rules <- c(zzz[(i_fam + 1L):(i_cur - 1L)], zzz[i_rpf:(i_end - 1L)])
curve_rules <- zzz[i_cur:(i_rpf - 1L)]
cat(sprintf("(a) curve: R/curve*.R %s; its compat rules in zzz.R %s\n",
            paste(loc(curve), collapse = "/"),
            paste(loc(curve_rules), collapse = "/")))
cat(sprintf("(b) RP: royston-parmar.R, rp-basis.R, rp-check.R %s; its compat rules %s\n",
            paste(loc(rp), collapse = "/"),
            paste(loc(rp_rules), collapse = "/")))
fam <- show("R/families.R")
a <- grep("^# --- Cox proportional hazards", fam)
b <- grep("^#' Call a family's log-CDF, passing", fam)
c1 <- grep("^#' The fitted baseline-hazard simplex of a", fam)
c2 <- grep("^#' The families whose response is an `n x K` matrix", fam)
d1 <- grep("^#' @section Cox proportional hazards:", fam)
d2 <- grep("^#' @section Quantile regression inference:", fam)
e1 <- grep("^#' @param df For `cox[(][)]`", fam)
e2 <- grep("^cox <- function", fam)
e3 <- e2 + which(fam[(e2 + 1L):length(fam)] == "}")[1L]
cox_fam <- c(fam[a:(b - 1L)], fam[c1:(c2 - 1L)], fam[d1:(d2 - 1L)],
             fam[(e1 - 1L):e3])
cat(sprintf("(c) cox(): families.R family, bases, cox_baseline(), docs, constructor %s\n",
            paste(loc(cox_fam), collapse = "/")))
comp <- show("R/compat.R")
ci <- grep("^  r[(]\"[^\"]+\", \"cox\",", comp)
cat(sprintf("    compat.R: %d hand rules naming cox, %d lines\n", length(ci),
            2L * length(ci)))
# the cox() branches elsewhere in core, each a block that read the
# baseline's own addition terms (Zbhaz, bhaz_chart, sbhaz_raw)
cnt <- function(path, pat) {
  x <- show(path)
  sum(grepl(pat, x))
}
for (p in c("R/objective.R", "R/fit.R", "R/priors.R", "R/frame.R",
            "R/predict.R", "R/sampling-api.R", "R/links-brms.R",
            "R/links.R", "R/compat.R", "R/families.R")) {
  cat(sprintf("    %-18s lines naming cox/bhaz: %d\n", p,
              cnt(p, "cox|bhaz")))
}

cat("\n== test time (sum of per-test real seconds, before run)\n")
O <- "dev/splitsurv-before-out"
rd <- function(pkg) {
  fs <- list.files(O, pattern = paste0("^", gsub(".", "[.]", pkg,
                                                  fixed = TRUE),
                                       "--.*[.]tsv$"), full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) {
    t <- utils::read.delim(f, stringsAsFactors = FALSE)
    if (!"real" %in% names(t)) return(NULL)
    t[, c("file", "test", "nb", "real")]
  }))
}
spl <- rd("frmtmb.spline")
core <- rd("frmtmb")
smp <- rd("frmtmb.sample")
tot_spl <- sum(spl$real)
tot_core <- sum(core$real)
pick <- function(t, files, tests = NULL) {
  s <- t$file %in% files
  if (!is.null(tests)) s <- s | (t$test %in% tests)
  t[s, ]
}
curve_t <- pick(spl, c("test-curve.R", "test-deriv.R", "test-difference.R",
                       "test-gratia.R", "test-illcond.R",
                       "test-new-levels.R", "test-span.R", "test-scale.R"))
curve_t <- rbind(curve_t, spl[spl$file == "test-surface.R" &
                                spl$test == "print methods say what was checked", ])
rp_t <- pick(spl, c("test-royston-parmar.R", "test-rp-floored.R",
                    "test-frailty.R"))
rp_t <- rbind(rp_t, spl[spl$file == "test-surface.R" &
                          spl$test != "print methods say what was checked", ])
cox_tests <- c(core$test[core$file == "test-famgaps.R" &
                           grepl("cox|M-spline|I-spline", core$test)])
cox_t <- rbind(core[core$file == "test-cox-chart.R", ],
               core[core$file == "test-famgaps.R" & core$test %in% cox_tests, ],
               core[core$file == "test-simulate-density.R" &
                      grepl("^cox[(][)] states why", core$test), ],
               core[core$file == "test-mv-gaps.R" &
                      grepl("extras are not per-response", core$test), ])
cat(sprintf("frmtmb.spline suite: %.1f s over %d tests\n", tot_spl,
            nrow(spl)))
cat(sprintf("(a) curve tests: %.1f s, %d tests, %.1f%% of the spline suite\n",
            sum(curve_t$real), nrow(curve_t),
            100 * sum(curve_t$real) / tot_spl))
cat(sprintf("(b) RP tests:    %.1f s, %d tests, %.1f%% of the spline suite\n",
            sum(rp_t$real), nrow(rp_t), 100 * sum(rp_t$real) / tot_spl))
cat(sprintf("core suite: %.1f s over %d tests\n", tot_core, nrow(core)))
cat(sprintf("(c) cox tests in core: %.1f s, %d tests, %.2f%% of the core suite\n",
            sum(cox_t$real), nrow(cox_t), 100 * sum(cox_t$real) / tot_core))
cd <- smp[smp$file == "test-cox-draws.R", ]
cat(sprintf("    and frmtmb.sample's test-cox-draws.R: %.1f s, %d tests, %.2f%% of that suite\n",
            sum(cd$real), nrow(cd), 100 * sum(cd$real) / sum(smp$real)))
print(aggregate(real ~ file, rbind(curve_t, rp_t, cox_t), sum))
