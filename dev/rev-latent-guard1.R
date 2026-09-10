# Reviewer, item 2.3, claim 1, part one: does the jitter actually move
# the starting values, and does hmm_starts_relist() write every outer
# component back?
#
# WHY THIS FIRST. hmm_starts_relist() skips a component whose
# fit[["estimates"]] entry is absent or the wrong length, silently and
# by `next`. If that skip fired for every component, every refit would
# start at the incumbent's own estimates, every refit would return the
# incumbent's own optimum, and the summary would report "1 distinct
# optimum, spread 0" while having tested nothing. That is a guard
# failing open, and it would look exactly like a reassuring result.
#
#   Rscript dev/rev-latent-guard1.R

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})

small_data <- function(seed = 4501L, N = 10L, Tl = 20L) {
  set.seed(seed)
  G <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, 2, byrow = TRUE)
  do.call(rbind, lapply(seq_len(N), function(id) {
    s <- integer(Tl); s[1L] <- 1L
    for (t in seq_len(Tl)[-1L]) s[t] <- sample.int(2L, 1L, prob = G[s[t - 1L], ])
    data.frame(id = id, t = seq_len(Tl),
               y = stats::rnorm(Tl, c(0, 3)[s], 0.6))
  }))
}

dd <- small_data()
fit <- frm(bf(y ~ 1), family = hmm(K = 2, gaussian(), time = t, group = id),
           data = dd)

par <- fit[["opt"]][["par"]]
cat("outer par names :", paste(names(par), collapse = " | "), "\n")
cat("estimates names :", paste(names(fit[["estimates"]]), collapse = " | "),
    "\n")
cat("lengths of estimates components:\n")
print(vapply(fit[["estimates"]], length, 1L))
cat("counts by outer par component:\n")
print(table(names(par)))

sc <- frmtmb.latent:::hmm_starts_scale(fit)
cat("\nscale_how:", sc$how, "\n")
cat("scale s   :", paste(format(sc$s, digits = 4), collapse = " "), "\n")
cat("all finite and positive:", all(is.finite(sc$s) & sc$s > 0), "\n")

## is confint()'s `est` column ANYTHING other than fit$opt$par?
ci <- suppressWarnings(stats::confint(fit))
cat("\nconfint() est identical to opt$par:",
    identical(as.numeric(ci[, "est"]), as.numeric(par)), "\n")
cat("max |est - opt$par|:",
    format(max(abs(as.numeric(ci[, "est"]) - as.numeric(par))),
           digits = 3), "\n")
cat("NOTE: confint.frmtmb_fit() builds its `est` column AS",
    "object$opt$par (R/confint.R:383),\n  so the alignment check in",
    "hmm_starts_scale() compares a vector with itself.\n")

## does the relist write every component?
set.seed(999)
v <- as.numeric(par) + stats::rnorm(length(par), 0, 2 * sc$s)
st <- frmtmb.latent:::hmm_starts_relist(fit, v)
cat("\nrelist: components in the returned start:",
    paste(names(st), collapse = " | "), "\n")
for (cp in unique(names(par))) {
  idx <- which(names(par) == cp)
  same <- isTRUE(all.equal(as.numeric(st[[cp]]),
                           as.numeric(fit[["estimates"]][[cp]])))
  cat(sprintf("  %-8s n_outer=%-3d n_est=%-3d moved=%s  max|move|=%s\n",
              cp, length(idx),
              length(fit[["estimates"]][[cp]]),
              if (same) "NO (SKIPPED)" else "yes",
              format(max(abs(as.numeric(st[[cp]]) -
                               as.numeric(fit[["estimates"]][[cp]]))),
                     digits = 3)))
}
moved <- vapply(unique(names(par)), function(cp) {
  !isTRUE(all.equal(as.numeric(st[[cp]]),
                    as.numeric(fit[["estimates"]][[cp]])))
}, TRUE)
cat("components actually perturbed:", sum(moved), "of", length(moved), "\n")

## and does a random-effect model behave the same?
cat("\n---- the d4 random-effect fit ----\n")
d4_data <- function() {
  set.seed(2026)
  K <- 2L; N <- 25L; Tg <- 30L
  G <- matrix(c(0.85, 0.15, 0.20, 0.80), 2, 2, byrow = TRUE)
  mu <- c(0, 3); sg <- c(0.6, 0.6); sd_b <- c(0.7, 0.5)
  stat <- local({
    A <- rbind(t(diag(K) - G), 1)
    drop(qr.solve(A, c(rep(0, K), 1)))
  })
  b <- cbind(stats::rnorm(N, 0, sd_b[1L]), stats::rnorm(N, 0, sd_b[2L]))
  d <- do.call(rbind, lapply(seq_len(N), function(g) {
    s <- integer(Tg); s[1L] <- sample.int(K, 1L, prob = stat)
    for (t in seq_len(Tg - 1L)) s[t + 1L] <- sample.int(K, 1L, prob = G[s[t], ])
    data.frame(ID = g, t = seq_len(Tg),
               y = stats::rnorm(Tg, mu[s] + b[g, s], sg[s]))
  }))
  d$gf <- factor(d$ID); d
}
d4 <- d4_data()
f4 <- frm(bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf)),
          family = hmm(K = 2, gaussian(), time = t, group = ID,
                       init = "stationary"),
          data = d4)
p4 <- f4[["opt"]][["par"]]
cat("outer par names :", paste(names(p4), collapse = " | "), "\n")
cat("estimates names :", paste(names(f4[["estimates"]]), collapse = " | "),
    "\n")
print(vapply(f4[["estimates"]], length, 1L))
sc4 <- frmtmb.latent:::hmm_starts_scale(f4)
cat("scale_how:", sc4$how, " s:",
    paste(format(sc4$s, digits = 4), collapse = " "), "\n")
set.seed(4101)
v4 <- as.numeric(p4) + stats::rnorm(length(p4), 0, 2 * sc4$s)
st4 <- frmtmb.latent:::hmm_starts_relist(f4, v4)
for (cp in unique(names(p4))) {
  idx <- which(names(p4) == cp)
  same <- isTRUE(all.equal(as.numeric(st4[[cp]]),
                           as.numeric(f4[["estimates"]][[cp]])))
  cat(sprintf("  %-8s n_outer=%-3d n_est=%-3d moved=%s\n", cp, length(idx),
              length(f4[["estimates"]][[cp]]),
              if (same) "NO (SKIPPED)" else "yes"))
}
cat("random-effect component `b` carried over unchanged (by design):",
    isTRUE(all.equal(as.numeric(st4[["b"]]),
                     as.numeric(f4[["estimates"]][["b"]]))), "\n")
