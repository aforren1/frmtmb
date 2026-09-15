## Re-check, items 2 and 3.
##
## 2. Is the null arm's second `V` an INDEPENDENT CONFIRMATION, or a
##    looser agreement being read alike? And was `V` computed the same
##    way there? A reference rung that is itself offset would move both
##    of that arm's residuals by the same amount, which is testable.
##
## 3. Is "the difference is the `id` rung's residual with the sign of
##    the comparison" an identity or a coincidence?
.libPaths(c("C:/Users/adf44/source/r/cohrev-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
source("dev/cohrev-read.R")
source("dev/coh-sim.R")

wide <- function(path) {
  d <- read_kv(path)
  keep <- names(which(table(d$seed) == length(unique(d$rung))))
  d <- d[as.character(d$seed) %in% keep, ]
  w <- reshape(d[, c("seed", "rung", "est")], idvar = "seed",
               timevar = "rung", direction = "wide")
  names(w) <- sub("^est[.]", "", names(w))
  w
}
w <- wide("dev/coh-recovery-main.tsv")
wn <- wide("dev/coh-recovery-null.tsv")

tr <- coupling_truth
fr <- seq_len(tr$n_freq) / tr$n_freq
v_bump <- mean((coupling_fbump(fr) - mean(coupling_fbump(fr)))^2)
att <- function(V) 1 / sqrt(1 + 0.346 * V)
pred <- function(V, V0 = 0) tr$b_cond * (att(V) - att(V0))

d_of <- function(tab, a, b) {
  x <- tab[[a]] - tab[[b]]
  c(mean = mean(x), se = sd(x) / sqrt(length(x)), n = length(x))
}
line <- function(lab, obs, p) {
  cat(sprintf("%-26s obs %+.5f se %.5f | pred %+.5f | resid %+.5f = %.2f se\n",
              lab, obs[["mean"]], obs[["se"]], p, obs[["mean"]] - p,
              (obs[["mean"]] - p) / obs[["se"]]))
}

cat("==== 2. the two arms side by side ====\n")
cat("main arm, reference `full`, V(full) = 0, n =", nrow(w), "\n")
mc <- d_of(w, "cond", "full")
ms <- d_of(w, "smooth", "full")
mi <- d_of(w, "id", "full")
line("cond", mc, pred(v_bump + tr$sd_id^2 + tr$sd_idcond^2))
line("smooth", ms, pred(tr$sd_id^2 + tr$sd_idcond^2))
line("id", mi, pred(tr$sd_idcond^2))
cat("null arm, reference `id`, V(id) = 0, n =", nrow(wn), "\n")
nc <- d_of(wn, "cond", "id")
ns <- d_of(wn, "smooth", "id")
line("cond", nc, pred(v_bump + tr$sd_id^2))
line("smooth", ns, pred(tr$sd_id^2))

cat("\nresiduals in ABSOLUTE terms, which is what 'confirmation' means here:\n")
cat(sprintf("  main arm: %+.5f and %+.5f, se %.5f\n",
            mc[["mean"]] - pred(v_bump + tr$sd_id^2 + tr$sd_idcond^2),
            ms[["mean"]] - pred(tr$sd_id^2 + tr$sd_idcond^2),
            ms[["se"]]))
cat(sprintf("  null arm: %+.5f and %+.5f, se %.5f\n",
            nc[["mean"]] - pred(v_bump + tr$sd_id^2),
            ns[["mean"]] - pred(tr$sd_id^2), ns[["se"]]))
cat(sprintf("  the null arm's se is %.1f times the main arm's, but its residual is %.1f times as large\n",
            ns[["se"]] / ms[["se"]],
            abs(ns[["mean"]] - pred(tr$sd_id^2)) /
              abs(ms[["mean"]] - pred(tr$sd_id^2 + tr$sd_idcond^2))))

## Both null-arm residuals are positive and about equal, which is the
## signature of an offset REFERENCE rather than a wrong prediction. In
## the null arm `full` is also correct, merely over-parameterized, so
## `full` minus `id` should be about zero if the reference is clean.
nf <- d_of(wn, "full", "id")
cat(sprintf("\nnull arm, `full` minus `id` (both correct there): %+.5f se %.5f, t %.2f\n",
            nf[["mean"]], nf[["se"]], nf[["mean"]] / nf[["se"]]))
cat(sprintf("if the reference carries %+.5f, the two residuals become %+.5f and %+.5f\n",
            nf[["mean"]],
            nc[["mean"]] - nf[["mean"]] - pred(v_bump + tr$sd_id^2),
            ns[["mean"]] - nf[["mean"]] - pred(tr$sd_id^2)))

cat("\n==== 3. identity or coincidence? ====\n")
p_smooth <- pred(tr$sd_id^2 + tr$sd_idcond^2)
p_id <- pred(tr$sd_idcond^2)
msi <- d_of(w, "smooth", "id")
r_smooth <- ms[["mean"]] - p_smooth
r_id <- mi[["mean"]] - p_id
r_si <- msi[["mean"]] - (p_smooth - p_id)
cat(sprintf("dropping (1 | id): pred %+.5f, obs %+.5f, resid %+.5f\n",
            p_smooth - p_id, msi[["mean"]], r_si))
cat(sprintf("the id rung's own residual: %+.5f\n", r_id))
cat(sprintf("the smooth rung's own residual: %+.5f\n", r_smooth))
cat(sprintf("r(smooth - id) + r(id) - r(smooth) = %.3g\n",
            r_si + r_id - r_smooth))
cat(sprintf("identical to machine precision: %s\n",
            isTRUE(all.equal(r_si, r_smooth - r_id, tolerance = 1e-12))))
cat("Both sides are differences of the SAME per-seed estimates, and\n")
cat("both the paired means and the predictions subtract, so the\n")
cat("residual subtracts too. It is arithmetic, not agreement: the\n")
cat("0.0002 the lane quotes IS the smooth rung's own residual.\n")
