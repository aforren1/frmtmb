# The versions of the seven extensions besides frmtmb.sample at 0.78.0,
# sourced by dev/rel078-resolve.R (local = TRUE). Decided by measurement
# (dev/round-20261029.md, "Versions and floors"): no file of spatial,
# survival, latent, eam, coupling, ode or learn changed, and every test
# file of theirs passes on the merged core (rellib-r16) with the count it
# has on the 0.77.0 core (rellib-r15), file by file, in the gated tier
# (dev/rel078-floor.sh, r078-floor-others, against r078-ref), and with
# 0.77.0's package totals in the plain tier: no bump, no new floor.
# So no DESCRIPTION but core's and frmtmb.sample's changes; checked here.
for (p in c("frmtmb.spatial", "frmtmb.survival", "frmtmb.latent",
            "frmtmb.eam", "frmtmb.coupling", "frmtmb.ode", "frmtmb.learn")) {
  f <- file.path("extensions", p, "DESCRIPTION")
  stopifnot(!f %in% c(merged, copied))
  h <- sub("\r$", "", system2("git", c("-C", git_dir, "show",
                                       paste0("HEAD:", f)), stdout = TRUE))
  stopifnot(identical(h, rd(f)))
}
cat("dev/rel078-versions.R: the seven other extensions keep their versions\n")
