# Reviewer, claim 6: record every ported brms families row (held,
# vacuous, message, what was caught) on one arm, and print the refusal
# messages of the rows that used to hold only on "could not find
# function".
#   Rscript dev/fams2-rev-port.R base|lane frmtmb|frmtmb.sample
args <- commandArgs(TRUE)
arm <- args[1]; pkg <- args[2]
base_lib <- c("C:/Users/adf44/source/r/rellib-r3",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lane") c("C:/Users/adf44/source/r/wt-fams2-lib",
                               base_lib) else base_lib)
sp <- paste0("C:/Users/adf44/AppData/Local/Temp/1/claude/",
             "c--Users-adf44-source-r-frmtmb/",
             "66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad")
root <- "C:/Users/adf44/source/r/frmtmb-wt-fams2"
tdir <- if (arm == "lane") {
  if (pkg == "frmtmb") file.path(root, "tests/testthat") else
    file.path(root, "extensions", pkg, "tests/testthat")
} else file.path(sp, "base", if (pkg == "frmtmb") "core" else "sample")
rec <- file.path(root, "dev/fams2-rev-out",
                 paste0("port-", arm, "-", pkg, ".tsv"))
unlink(rec)
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_RECORD = rec, FRMTMB_BRMSPORT_PKG = pkg,
           FRMTMB_STAN_CACHE = file.path(sp, "stan-cache"))
suppressPackageStartupMessages(library(testthat))
invisible(test_file(file.path(tdir, "test-brms-suite-families.R"),
                    package = pkg, env = test_env(pkg), reporter = "silent",
                    load_package = "installed"))
r <- utils::read.delim(rec, header = FALSE, quote = "",
                       col.names = c("kind", "pkg", "id", "verdict", "held",
                                     "vacuous", "msg", "caught", "raw"))
cat(arm, pkg, ": rows recorded", nrow(r), "\n")
flip <- r[r$verdict != "pass" & r$held == "TRUE", ]
cat("verdict not pass but HELD:", nrow(flip), "\n")
cat(paste(flip$id, collapse = " "), "\n")
stay <- r[r$id %in% paste0("families:", c(107, 108, 109, 117, 118, 122)), ]
print(stay[, c("id", "verdict", "held", "vacuous", "msg")], right = FALSE)
if (arm == "lane") {
  library(pkg, character.only = TRUE)
  cat("\nThe refusals rows 57-59, 65-67, 71, 74 now hold on:\n")
  calls <- list(
    `57` = quote(zero_inflated_beta_binomial('sqrt')),
    `58` = quote(zero_inflated_beta_binomial(link_phi = 'logit')),
    `59` = quote(zero_inflated_beta_binomial(link_zi = 'log')),
    `65` = quote(hurdle_cumulative(link = "log")$link),
    `66` = quote(hurdle_cumulative(link_hu = "probit")$link_hu),
    `67` = quote(hurdle_cumulative(link_disc = "logit")$link_disc),
    `71` = quote(xbeta("1/mu")),
    `74` = quote(xbeta(link_phi = "sqrt")$link_phi))
  for (nm in names(calls)) {
    m <- tryCatch({ eval(calls[[nm]]); "NO ERROR" },
                  error = function(e) paste0("[", class(e)[1], "] ",
                                             conditionMessage(e)))
    cat(sprintf("  %s: %s\n", nm, substr(m, 1, 150)))
  }
  cat("\nAbsent-condition twins (brms accepts these):\n")
  for (cl in list(quote(zero_inflated_beta_binomial('probit')$link),
                  quote(zero_inflated_beta_binomial(link_zi = 'identity')$link_zi),
                  quote(hurdle_cumulative(link_hu = "identity")$link_hu),
                  quote(hurdle_cumulative(link_disc = "identity")$link_disc),
                  quote(hurdle_cumulative("cloglog")$link),
                  quote(xbeta(link_phi = "identity")$link_phi))) {
    m <- tryCatch(format(eval(cl)), error = function(e) paste("ERROR", conditionMessage(e)))
    cat(sprintf("  %s -> %s\n", deparse1(cl), substr(m, 1, 100)))
  }
  cat("\nbrms 2.23.0 on the same calls:\n")
  for (cl in list(quote(brms::hurdle_cumulative(link_hu = "identity")$link_hu),
                  quote(brms::hurdle_cumulative(link_disc = "identity")$link_disc),
                  quote(brms::zero_inflated_beta_binomial(link_zi = 'identity')$link_zi),
                  quote(brms::xbeta(link_phi = "identity")$link_phi))) {
    m <- tryCatch(format(eval(cl)), error = function(e) paste("ERROR", conditionMessage(e)))
    cat(sprintf("  %s -> %s\n", deparse1(cl), substr(m, 1, 100)))
  }
}
