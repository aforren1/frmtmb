# The draws stand-in for brms:::brmsfit_example1..6, used by the
# frmtmb.sample half of the brms suite port (item 2.6b). It lives here
# and not in helper-brms-suite.R because core does not suggest
# frmtmb.sample, and dev/brmsport-gen.R copies that helper verbatim.
# brms's example fits hold 25 post-warmup draws of one chain of 75
# iterations, and the suite asserts those counts (ndraws 25, nsamples
# with warmup 75), so the draws stand-in is drawn to the same design.
# The seed only fixes the draws; no assertion reads their values.
brms_fixture_draws <- function(k) {
  key <- paste0("draws", k)
  hit <- brms_port_state$fixtures[[key]]
  if (!is.null(hit)) return(hit)
  fit <- brms_fixture(k)
  ds <- withCallingHandlers(
    frmtmb.sample::frm_sample(fit, chains = 1, iter = 75, warmup = 50,
                              seed = 20260917, refresh = 0),
    warning = function(w) invokeRestart("muffleWarning"),
    message = function(m) invokeRestart("muffleMessage")
  )
  brms_port_state$fixtures[[key]] <- ds
  ds
}
