# Regenerate inst/extdata/habit-hierarchical.csv, the precomputed group-level
# model comparison quoted in vignette("habit").
#
# Six fits, one habit and one no-habit model per group. They are precomputed
# because together they take a few minutes, which is longer than a vignette
# should spend. Each is run from several starting values: the likelihood
# surface for the habit model has more than one local optimum, and the
# minimal-practice group in particular lands in different basins depending on
# where it starts.
#
# Run from the package root:
#   Rscript data-raw/habit_hierarchical.R

library(frmtmb)
data(habit_prep, package = "frmtmb")

dat <- subset(habit_prep, remapped & !is.na(response))
dat$response <- factor(dat$response, levels = c("correct", "habit", "other"))

habit_hier <-
  bf(response ~ 1, family = categorical()) +
  nlf(PhiA ~ pnorm((prep_time - muA) / sgA)) +
  nlf(PhiB ~ pnorm((prep_time - muB) / sgB)) +
  nlf(sgA ~ exp(lsgA)) + nlf(sgB ~ exp(lsgB)) +
  nlf(qB ~ 0.5 + 0.4999 / (1 + exp(-lqB))) +
  nlf(qI ~ 0.499 / (1 + exp(-lqI))) +
  nlf(c1 ~ qI*(1-PhiA)*(1-PhiB) + ((1-0.99)/3)*PhiA*(1-PhiB) + qB*PhiB) +
  nlf(c2 ~ qI*(1-PhiA)*(1-PhiB) + 0.99*PhiA*(1-PhiB) + ((1-qB)/3)*PhiB) +
  nlf(c3 ~ (0.5-qI)*(1-PhiA)*(1-PhiB) + ((1-0.99)/3)*PhiA*(1-PhiB) +
           ((1-qB)/3)*PhiB) +
  nlf(muhabit ~ log(c2 / c1)) + nlf(muother ~ log(2 * c3 / c1)) +
  lf(muA ~ 1 + (1 | p | participant), muB ~ 1 + (1 | p | participant),
     lsgA ~ 1 + (1 | p | participant), lsgB ~ 1 + (1 | p | participant),
     lqB ~ 1, lqI ~ 1)

nohabit_hier <-
  bf(response ~ 1, family = categorical()) +
  nlf(PhiB ~ pnorm((prep_time - muB) / sgB)) +
  nlf(sgB ~ exp(lsgB)) +
  nlf(qB ~ 0.5 + 0.4999 / (1 + exp(-lqB))) +
  nlf(qI ~ 0.499 / (1 + exp(-lqI))) +
  nlf(c1 ~ qI*(1-PhiB) + qB*PhiB) +
  nlf(c2 ~ qI*(1-PhiB) + ((1-qB)/3)*PhiB) +
  nlf(c3 ~ (0.5-qI)*(1-PhiB) + ((1-qB)/3)*PhiB) +
  nlf(muhabit ~ log(c2 / c1)) + nlf(muother ~ log(2 * c3 / c1)) +
  lf(muB ~ 1 + (1 | p | participant), lsgB ~ 1 + (1 | p | participant),
     lqB ~ 1, lqI ~ 1)

nm_h <- c("muA_(Intercept)", "lsgA_(Intercept)", "muB_(Intercept)",
          "lsgB_(Intercept)", "lqB_(Intercept)", "lqI_(Intercept)")
nm_n <- c("muB_(Intercept)", "lsgB_(Intercept)", "lqB_(Intercept)",
          "lqI_(Intercept)")
starts_h <- list(c(0.4, log(0.08), 0.5, log(0.08), 2, 0),
                 c(0.6, log(0.15), 0.45, log(0.05), 1.5, 0.5),
                 c(1.2, log(0.30), 0.5, log(0.08), 2, 0))
starts_n <- list(c(0.5, log(0.08), 2, 0), c(0.45, log(0.05), 1.5, 0.5),
                 c(0.55, log(0.15), 2.5, -0.3))

best_of <- function(form, d, starts, nms) {
  best <- NULL; best_ll <- -Inf; n_ok <- 0
  for (v in starts) {
    f <- try(suppressWarnings(
      frm(form, data = d, start = list(beta = setNames(v, nms)))), silent = TRUE)
    if (inherits(f, "try-error")) next
    n_ok <- n_ok + 1
    ll <- as.numeric(logLik(f))
    if (is.finite(ll) && ll > best_ll) { best_ll <- ll; best <- f }
  }
  list(fit = best, n_ok = n_ok)
}

rows <- lapply(c("minimal", "4day", "20day"), function(g) {
  d <- droplevels(dat[dat$group == g, ])
  h <- best_of(habit_hier, d, starts_h, nm_h)
  n <- best_of(nohabit_hier, d, starts_n, nm_n)
  llh <- logLik(h$fit); lln <- logLik(n$fit)
  data.frame(group = g, n_trials = nrow(d), n_subj = nlevels(d$participant),
             logLik_habit = as.numeric(llh), df_habit = attr(llh, "df"),
             AIC_habit = AIC(h$fit),
             logLik_nohabit = as.numeric(lln), df_nohabit = attr(lln, "df"),
             AIC_nohabit = AIC(n$fit),
             dAIC = AIC(n$fit) - AIC(h$fit),
             starts_ok_habit = h$n_ok, starts_ok_nohabit = n$n_ok)
})

write.csv(do.call(rbind, rows), "inst/extdata/habit-hierarchical.csv",
          row.names = FALSE)
