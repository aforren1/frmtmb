# L1 of the frm_sample() route against brms's own program, ordinary noise.
#
#   Rscript dev/gpgap-sampler-l1.R "<spec key>"
#
# The sampler's objective, its tape non-centered at its own nugget (1e-12 on
# the correlation, frmtmb.sample's gp_sampler_nugget), at frm()'s estimate
# with z through its own factor; brms's program at z through that same kernel
# in brms's order. Priors are dropped: the bare likelihood and the block on
# both sides. Also the same reading at brms's line for frmtmb's factor
# (`frmline`), the identity check that nothing else moved: there the two
# programs differ by the kernel line only, so brms's rewritten line must read
# 0. Lane build.
key <- commandArgs(trailingOnly = TRUE)[1]
Sys.setenv(GPGAP_LIB = paste0("C:/Users/adf44/source/r/wt-gpgap-lib;",
                              "C:/Users/adf44/source/r/rellib-r16"))
source("dev/gpgap-lib.R")
suppressPackageStartupMessages(library(frmtmb.sample))
E$key <- key
local(envir = E, {
  spec <- eq_specs[[key]]()
  spec$prior <- NULL
  spec$l1_code <- NULL
  s <- gpgap_setup(spec)
  nug <- utils::getFromNamespace("gp_sampler_nugget", "frmtmb.sample")
  fit <- s$fit
  ib <- which(vapply(fit$frame$re_blocks, function(b) {
    identical(b$covstruct, "gp")
  }, NA))
  # the sampler's own block: its nugget and its factor's order
  fs <- utils::getFromNamespace("gp_sample_nugget", "frmtmb.sample")(fit, ib)
  fr <- fs$frame
  fr$ncp_blocks <- ib
  stopifnot(identical(fr$re_blocks[[ib]]$gp_nugget, nug))
  fn <- frmtmb::build_objective(fr)
  bk <- fr$re_blocks[[ib]]
  th <- fit$estimates$theta[bk$theta_idx]
  L <- frmtmb:::ncp_block_chol(bk, th)
  pz <- fit$estimates
  pz$b[bk$b_idx] <- as.numeric(solve(L, fit$estimates$b[bk$b_idx]))
  # a point of the sampler's own typical set, not frm()'s field through the
  # smaller nugget, whose z would be the size of the near-null directions: z
  # from frm()'s factor in brms's order, and the field that z makes through
  # the sampler's factor in that order. The tape is sorted, so its z is read
  # back through its own factor, and the likelihood is the tape less that
  # z's standard normal.
  io <- brms_gp_order(s$sdat, s$zn)
  zb <- as.numeric(solve(t(chol(s$kb)), s$fb))
  Ls <- t(chol(tcrossprod(L)[io, io]))
  f <- numeric(length(zb))
  f[io] <- as.numeric(Ls %*% zb)
  pz$b[bk$b_idx] <- as.numeric(solve(L, f))
  ll_tape <- -fn(pz) - sum(stats::dnorm(pz$b[bk$b_idx], log = TRUE))
  ours <- ll_tape + sum(stats::dnorm(zb, log = TRUE))
  p <- s$pars
  p[[s$zn]] <- array(zb, length(zb))
  sf0 <- gpgap_sf(s$code, s$sdat)
  own <- gpgap_lp(sf0, p) - ours
  # the same tape at frmtmb's 1e-6, against brms with that line rewritten
  fr2 <- fit$frame
  fr2$ncp_blocks <- ib
  fn2 <- frmtmb::build_objective(fr2)
  bk2 <- fr2$re_blocks[[ib]]
  L2 <- frmtmb:::ncp_block_chol(bk2, th)
  pz2 <- fit$estimates
  pz2$b[bk2$b_idx] <- as.numeric(solve(L2, fit$estimates$b[bk2$b_idx]))
  Lb2 <- t(chol(s$kb))
  p2 <- s$pars
  p2[[s$zn]] <- array(as.numeric(solve(Lb2, s$fb)), length(s$fb))
  sf1 <- gpgap_sf(gpgap_frm_nugget(s$code), s$sdat)
  check <- gpgap_lp(sf1, p2) - (-fn2(pz2))
  cat(sprintf(paste0("GPSAMPL1 %s: sdgp %.4g, ours %.3f, tol %.2g | sampler ",
                     "tape against brms's own line %.3e | tape at 1e-6 against ",
                     "brms rewritten %.3e\n"),
              key, s$sdgp, ours, 1e-6 * abs(ours), own, check))
})
