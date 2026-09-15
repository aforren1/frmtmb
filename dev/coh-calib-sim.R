## The small design the shipped assertion uses, and the model-free
## reference standard error it is measured against. Sourced by
## dev/coh-calib2.R and copied into the package's test file once the
## numbers below settled.

## No frequency axis: the contrast is what is being measured, so the
## rows inside a subject-by-condition cell are replicate draws and
## nothing in the truth depends on which one it is.
coh_flat <- function(seed, n_sub, n_rep, sd_idcond, sd_id = 0.35,
                     b0 = -0.6, b_cond = 0.5, nseg = 8L) {
  set.seed(seed)
  d <- expand.grid(rep = seq_len(n_rep), cond = factor(c("a", "b")),
                   id = factor(seq_len(n_sub)))
  u_id <- stats::rnorm(n_sub, 0, sd_id)
  u_ic <- stats::rnorm(n_sub * 2L, 0, sd_idcond)
  ic <- as.integer(d$id) + n_sub * (as.integer(d$cond) - 1L)
  eta <- b0 + b_cond * (d$cond == "b") +
    u_id[as.integer(d$id)] + u_ic[ic]
  w <- coupling_draw(rep(exp(0.3), nrow(d)), rep(exp(0.1), nrow(d)),
                     stats::plogis(eta), rep(0.4, nrow(d)), nseg)
  cbind(d, w)
}

## The two-stage standard error: pool each cell, take its coherence on
## the logit scale, difference the conditions within a subject, and read
## the standard error off the between-subject spread of the
## differences. It carries the subject-by-condition spread because that
## spread IS the spread of the differences.
coh_two_stage <- function(d) {
  cell <- function(i, cd) {
    k <- d$id == i & d$cond == cd
    w11 <- sum(d$w11[k]); w22 <- sum(d$w22[k])
    r <- sum(d$w12r[k]); im <- sum(d$w12i[k])
    stats::qlogis((r^2 + im^2) / (w11 * w22))
  }
  ids <- levels(d$id)
  dif <- vapply(ids, function(i) cell(i, "b") - cell(i, "a"),
                numeric(1))
  list(est = mean(dif), se = stats::sd(dif) / sqrt(length(dif)),
       n = length(dif))
}
