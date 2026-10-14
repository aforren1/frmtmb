# Reviewer round 2: distributions:238 is recorded as 8 divergences.
# Which of its 8 assertion lines hold for frmtmb over EVERY iteration of
# brms's loops? Seed 235 (the generated test's local_seed), the lane's
# adapter brms_adapt_dord(), brms's own reference helpers.
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
`%||%` <- function(a, b) if (is.null(a)) b else a
source("dev/brmsport-bin2-adapters.R")
sd <- "dev/brms-suite/brms/tests/testthat/helpers"
source(file.path(sd, "inv_link_ordinal_ch.R"))
source(file.path(sd, "simopts_catlike_oneobs.R"))
eqv <- function(a, b) isTRUE(all.equal(unname(a), unname(b),
                                       check.attributes = FALSE))
tally <- list()
note <- function(id, ok) tally[[id]] <<- c(tally[[id]], ok)
set.seed(235)
suppressWarnings(
for (ndraws in ndraws_vec) for (ncat in ncat_vec) {
  thres_test <- matrix(rnorm(ndraws * (ncat - 1)), nrow = ndraws)
  eta_test_list <- list(rnorm(ndraws),
                        matrix(rnorm(ndraws * (ncat - 1)), nrow = ndraws))
  for (eta_test in eta_test_list) {
    thres_eta <- if (is.matrix(eta_test)) thres_test - eta_test else
      sweep(thres_test, 1, as.array(eta_test))
    eta_thres <- if (is.matrix(eta_test)) eta_test - thres_test else
      sweep(-thres_test, 1, as.array(eta_test), FUN = "+")
    for (link in c("logit", "probit", "cauchit", "cloglog")) {
      chk <- list(cumulative = inv_link_cumulative_ch(thres_eta, link),
                  sratio = inv_link_sratio_ch(thres_eta, link),
                  cratio = inv_link_cratio_ch(eta_thres, link),
                  acat = inv_link_acat_ch(eta_thres, link))
      for (f in names(chk)) {
        d <- brms_adapt_dord(f, seq_len(ncat), eta_test, thres_test,
                             link = link)
        note(paste(f, "equivalent"), eqv(d, chk[[f]]))
        note(paste(f, "dim"), identical(dim(d), c(as.integer(ndraws),
                                                  as.integer(ncat))) ||
               isTRUE(all.equal(dim(d), c(ndraws, ncat))))
      }
    }
  }
})
for (id in names(tally)) {
  cat(sprintf("%-24s holds in %d of %d iterations\n", id,
              sum(tally[[id]]), length(tally[[id]])))
}
