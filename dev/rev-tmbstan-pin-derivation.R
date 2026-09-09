# Reviewer demonstration for lane tmbstan, finding R1.
#
# The pin step in .github/workflows/check-frmtmb-sample.yaml derives its
# dated snapshot from `getOption("repos")[["CRAN"]]`. On a runner
# configured by r-lib/actions/setup-r@v2 with use-public-rspm: true the
# CRAN entry is NOT the RSPM URL. setup-r writes an .Rprofile holding
#
#     repos <- c(RSPM = <the packagemanager URL>, CRAN = <the mirror>)
#
# and the mirror defaults to https://cran.rstudio.com because the
# action's `cran` input has no default (setup-r/action.yml) and neither
# the CRAN nor the RSPM environment variable is set by this workflow.
# The runner prints the resulting table itself; run 34322... of
# check-frmtmb-sample, step "Run r-lib/actions/setup-r-dependencies@v2",
# "Repo status":
#
#     1  RSPM  https://packagemanager.posit.co/cran/__linux__/noble/latest
#     2  CRAN  https://cran.rstudio.com
#
# This script replays the step's own arithmetic on both entries.
say <- function(label, repo) {
  pin <- sub("/latest/?$", "/2026-09-01", repo)
  cat(sprintf("%-6s in : %s\n", label, repo))
  cat(sprintf("%-6s out: %s\n", label, pin))
  cat(sprintf("%-6s no-op (step calls stop()): %s\n\n", label,
              identical(pin, repo)))
}

say("CRAN", "https://cran.rstudio.com")
say("RSPM", "https://packagemanager.posit.co/cran/__linux__/noble/latest")
say("RSPM/", "https://packagemanager.posit.co/cran/__linux__/noble/latest/")

# and the exit status the runner sees for the no-op branch
repo <- "https://cran.rstudio.com"
pin <- sub("/latest/?$", "/2026-09-01", repo)
if (identical(pin, repo)) {
  cat("the step reaches its stop() branch, so the step exits non-zero\n")
  cat("and the job fails before install.packages() is ever called\n")
}
