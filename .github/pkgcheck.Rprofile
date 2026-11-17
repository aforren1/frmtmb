# Sourced by the pkgcheck-action container via R_PROFILE_USER (set in
# .github/workflows/pkgcheck.yaml). Every dependency of the core package
# is on CRAN; the ODE extension is the one r-universe consumer and it is
# not checked here.

# The container resolves dependencies with pak::lockfile_create(), which
# reads getOption("repos"). Left at CRAN's source repository it built
# all 222 packages from source, rstan and BH included, which took 20 of
# the job's 60 minutes. Posit Package Manager serves prebuilt binaries
# for this container's Ubuntu 24.04 (noble); a user profile runs after
# the site profile, so this setting wins.
options(repos = c(
  CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"
))

# belt and suspenders for the workflow env (see pkgcheck.yaml)
Sys.setenv(FRMTMB_SAMPLER_GATES = "false")

# Two defects of pkgcheck 0.3.4.5 fail this job whatever the package's
# state (dev/upstream-bugs.md, pkgcheck-1 and pkgcheck-2). They are
# patched here, in this job's R session only, until upstream fixes them;
# each patch does nothing once the defect is gone.
local({
  # pkgcheck-2: the CI check reads the workflow runs from the GitHub API
  # with no token and swallows any error, so a failed anonymous request
  # reports "Package has no continuous integration checks" although
  # R-CMD-check and test-coverage passed. The container gets the job's
  # GITHUB_TOKEN; send it with requests to the API.
  auth_github <- function(...) {
    tok <- Sys.getenv("GITHUB_TOKEN")
    if (!nzchar(tok)) return(invisible())
    request <- get("request", envir = asNamespace("httr2"))
    utils::assignInNamespace("request", function(base_url) {
      req <- request(base_url)
      if (startsWith(base_url, "https://api.github.com/")) {
        req <- httr2::req_auth_bearer_token(req, tok)
      }
      req
    }, ns = "httr2")
  }
  # pkgcheck-1 (pkgcheck issue #444): when the CI check passes,
  # checks_to_markdown() returns a list, and the action then stops at
  # writeLines() with "can only write character objects".
  flat_markdown <- function(...) {
    md <- get("checks_to_markdown", envir = asNamespace("pkgcheck"))
    utils::assignInNamespace("checks_to_markdown", function(...) {
      out <- md(...)
      if (is.list(out)) unlist(out, use.names = FALSE) else out
    }, ns = "pkgcheck")
  }
  if (isNamespaceLoaded("httr2")) auth_github() else
    setHook(packageEvent("httr2", "onLoad"), auth_github)
  if (isNamespaceLoaded("pkgcheck")) flat_markdown() else
    setHook(packageEvent("pkgcheck", "onLoad"), flat_markdown)
})
