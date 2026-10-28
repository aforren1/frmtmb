# The versions of the seven extensions besides frmtmb.sample at 0.77.0,
# sourced by dev/rel077-resolve.R (local = TRUE). Decided by measurement
# (dev/round-20261028.md, "Versions and floors"):
# - frmtmb.spatial: no R file changed; sigmainn changed its copy of the
#   test helper helper-brms.R (the 0.76.0 innovation-sd correction of
#   stan_pars_from_fit() went with the change of sigma). Its 8 test files
#   from this tree pass on rellib-r14 (frmtmb 0.76.0) with the counts
#   they have on rellib-r15 (dev/rel077-floor.sh), so its tests do not
#   need the new core: no bump, floor frmtmb (>= 0.76.0) stands.
# - survival, latent, eam, coupling, ode, learn: no file changed, and
#   every test file of theirs passes on the merged core in both tiers
#   with its 0.76.0 count: no bump.
# So no DESCRIPTION but core's and frmtmb.sample's changes; checked here.
for (p in c("frmtmb.spatial", "frmtmb.survival", "frmtmb.latent",
            "frmtmb.eam", "frmtmb.coupling", "frmtmb.ode", "frmtmb.learn")) {
  f <- file.path("extensions", p, "DESCRIPTION")
  stopifnot(!f %in% c(merged, copied))
  h <- sub("\r$", "", system2("git", c("-C", git_dir, "show",
                                       paste0("HEAD:", f)), stdout = TRUE))
  stopifnot(identical(h, rd(f)))
}
cat("dev/rel077-versions.R: the seven other extensions keep their versions\n")
