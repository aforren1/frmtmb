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
