# Sourced by the pkgcheck-action container via R_PROFILE_USER (set in
# .github/workflows/pkgcheck.yaml). The action installs the package's
# dependencies with pak, which cannot read Additional_repositories from
# DESCRIPTION; RTMBode (Suggests) only exists on r-universe. At
# user-profile time getOption("repos") can still be NULL (the CRAN
# default is applied later in startup), so a CRAN entry is set
# explicitly rather than relying on an append.
local({
  r <- getOption("repos")
  if (!length(r) || identical(unname(r), "@CRAN@")) {
    r <- c(CRAN = "https://cloud.r-project.org")
  }
  options(repos = c(r, kaskr = "https://kaskr.r-universe.dev"))
})
