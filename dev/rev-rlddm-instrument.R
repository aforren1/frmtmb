# REVIEW, item 1.0b, attack 2: the instrument the population
# non-decision time is measured with, on both scale rows.
#
#   Rscript dev/rev-rlddm-instrument.R
#
# No seed: arithmetic on the two rows' own recorded fields,
# dev/rlddm-scripts/rlddm-scale-after.tsv and
# dev/ndt-scripts/ndt-scale-r1.tsv.

f <- function(lab, ns, sdlog, mu, se, pop, sub_true, sub_hat) {
  e_ndt <- mu * exp(sdlog^2 / 2)
  sd_nat <- mu * sqrt(exp(sdlog^2) * (exp(sdlog^2) - 1))
  sem <- sd_nat / sqrt(ns)
  cat(sprintf("%-6s ns=%3d  E[ndt]=%.6f  sd_nat=%.6f\n", lab, ns,
              e_ndt, sd_nat))
  cat(sprintf("%-6s sem of the DRAW %.6f, se the row reports %.6f, ",
              "", sem, se))
  cat(sprintf("ratio %.2f\n", sem / se))
  cat(sprintf("%-6s draw %+.2f ms vs E, fit %+.2f ms, conversion %+.2f ms",
              "", 1000 * (sub_true - e_ndt), 1000 * (sub_hat - sub_true),
              1000 * (pop - sub_hat)))
  cat(sprintf(", z = %.2f\n\n", abs(pop - mu) / se))
  invisible(list(e = e_ndt, sem = sem))
}

a <- f("learn", 100, 0.15, 0.25, 0.00175853, 0.25929, 0.254744, 0.257234)
b <- f("eam", 30, 0.12, 0.25, 0.0019639, 0.246867, 0.243452, 0.246146)

# how much draw luck the eam row has before its own z < 4 breaks, with
# the fit and the conversion held at what this seed produced
off <- 0.246867 - 0.243452
brk <- 0.25 + 4 * 0.0019639
need <- brk - off
cat(sprintf("eam: z < 4 breaks once the realized draw mean passes %.6f\n",
            need))
cat(sprintf("     that is E[ndt] %+.2f ms, %+.2f sem of the draw\n",
            1000 * (need - b$e), (need - b$e) / b$sem))
cat(sprintf("     P(a draw mean at least that high) = %.3f\n",
            stats::pnorm(need, b$e, b$sem, lower.tail = FALSE)))

off2 <- 0.25929 - 0.254744
brk2 <- 0.25 + 4 * 0.00175853
need2 <- brk2 - off2
cat(sprintf("learn: z < 4 would break past %.6f, E[ndt] %+.2f ms, ",
            need2, 1000 * (need2 - a$e)))
cat(sprintf("%+.2f sem, P = %.3f\n", (need2 - a$e) / a$sem,
            stats::pnorm(need2, a$e, a$sem, lower.tail = FALSE)))
