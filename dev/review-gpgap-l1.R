# Short check of fix round 1, item 5: is the sampler's L1 against brms's own line
# "relative against absolute"? At sdgp = 1 the two nuggets are one number, so a
# gap that survives there is not the nugget. Three readings at each sdgp, one
# typical z: tape less brms's own line; tape less brms's program with its line
# rewritten to the sampler's relative 1e-12; own less rewritten, Stan on both
# sides, so the same arithmetic. Same point construction as the lane's
# dev/gpgap-sampler-l1.R.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-l1.R <log> "<spec key>"
key <- commandArgs(trailingOnly = TRUE)[1]
Sys.setenv(GPGAP_LIB = paste0("C:/Users/adf44/source/r/wt-gpgap-lib;",
                              "C:/Users/adf44/source/r/rellib-r16"))
cache <- Sys.getenv("FRMTMB_STAN_CACHE")
source("dev/gpgap-lib.R")
Sys.setenv(FRMTMB_STAN_CACHE = cache)
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
  fr <- fit$frame
  fr$ncp_blocks <- ib
  fr$re_blocks[[ib]]$gp_nugget <- nug
  fn <- frmtmb::build_objective(fr)
  bk <- fr$re_blocks[[ib]]
  th0 <- fit$estimates$theta[bk$theta_idx]
  io <- brms_gp_order(s$sdat, s$zn)
  # one typical z: frm()'s field through frm()'s factor in brms's order
  zb <- as.numeric(solve(t(chol(s$kb)), s$fb))
  sf_own <- gpgap_sf(s$code, s$sdat)
  rel <- gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-12 * square(sdgp);", s$code,
              fixed = TRUE)
  stopifnot(!identical(rel, s$code))
  sf_rel <- gpgap_sf(rel, s$sdat)
  sdn <- paste0("sdgp_", s$gsfx)
  read <- function(sd) {
    th <- th0
    th[[1]] <- log(sd)
    L <- frmtmb:::ncp_block_chol(bk, th)
    # the field the tape's factor makes of this z in brms's order
    Ls <- t(chol(tcrossprod(L)[io, io]))
    f <- numeric(length(zb))
    f[io] <- as.numeric(Ls %*% zb)
    pz <- fit$estimates
    pz$theta[bk$theta_idx] <- th
    pz$b[bk$b_idx] <- as.numeric(solve(L, f))
    ll_tape <- -fn(pz) - sum(stats::dnorm(pz$b[bk$b_idx], log = TRUE))
    ours <- ll_tape + sum(stats::dnorm(zb, log = TRUE))
    p <- s$pars
    p[[s$zn]] <- array(zb, length(zb))
    p[[sdn]] <- array(sd, 1L)
    a <- gpgap_lp(sf_own, p)
    b <- gpgap_lp(sf_rel, p)
    # the same z through two factors of one kernel, by R, on the tape's sort and
    # in brms's order: how far rounding moves the field at this condition
    K <- exp(2 * th[[1]]) * tcrossprod(L / exp(th[[1]]))
    Lb <- t(chol(K[io, io]))
    fb2 <- numeric(length(zb))
    fb2[io] <- as.numeric(Lb %*% zb)
    c(tape_own = a - ours, tape_rel = b - ours, own_rel = a - b,
      field = max(abs(f - fb2)), ours = ours)
  }
  sds <- c(s$sdgp, 1, 1 + 1e-6, 1 - 1e-6, 1.001, 0.999, 0.95, 1.05, 1.5)
  for (sd in sds) {
    v <- read(sd)
    cat(sprintf(paste0("L1S %s sdgp %.7g | tape - brms own %10.3e | tape - brms rel 1e-12 ",
                       "%10.3e | own - rel, Stan both %10.3e | ours %.3f\n"),
                key, sd, v[["tape_own"]], v[["tape_rel"]], v[["own_rel"]], v[["ours"]]))
  }
})
