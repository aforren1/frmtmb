# Lane `learnhier`: how does the nuisance term move with the trial
# count? NO FITS AT ALL.
#
#   Rscript dev/learnhier-mvsn.R [n_seeds]
#
# The fitted deviation under `ndt_group(id)` is exactly
# `log(ndt_i) - log(m_i)`, where `m_i` is that learner's smallest
# first-passage time. `sd(log ndt)` is fixed by the design at 0.15.
# Everything that moves is `m`, and `m` is a property of SIMULATED DATA,
# so this needs no model, no optimizer and no fit.
#
# The arms are PAIRED on their learners. lh_rlddm_data() draws the block
# from `set.seed(seed + 2)` before it simulates anything, so one seed
# gives the same hundred learners with the same true parameters at every
# trial count. Only the data they produce, and therefore their minima,
# differ.
#
# The committed direction is in dev/learnhier-findings.md under
# "Test 2's prediction REVERSED": more trials do NOT repair the
# parameterization, so the nuisance share does not fall. What is certain
# is only that the LEVEL of `m` falls, since a minimum over more draws
# cannot rise; the direction of `sd(log m)` is not obvious on paper
# because the Wiener first-passage density near zero has an essential
# singularity rather than a power-law tail.
source("dev/learnhier-env.R")
source("dev/learnhier-sim.R")
suppressPackageStartupMessages(library(frmtmb.learn))

a <- commandArgs(trailingOnly = TRUE)
NS <- if (length(a)) as.integer(a[[1L]]) else 12L
seeds <- 20261000L + 10L * seq_len(NS)
counts <- c(100L, 200L, 400L)

rows <- list()
for (nt in counts) {
  for (sd_i in seeds) {
    d <- lh_rlddm_data(sd_i, ns = 100L, nt = nt)
    nd <- as.numeric(attr(d, "ndt_subject"))
    fl <- as.numeric(attr(d, "own_floor"))
    m <- fl - nd
    eta <- log(nd) - log(m)
    dv <- attr(d, "dev_drawn")
    rows[[length(rows) + 1L]] <- data.frame(
      nt = nt, seed = sd_i,
      mean_log_m = mean(log(m)),
      sd_log_m = stats::sd(log(m)),
      sd_log_ndt = stats::sd(log(nd)),
      sd_eta = stats::sd(eta),
      share_min = stats::sd(log(m))^2 / stats::sd(eta)^2,
      cor_eta_bs = stats::cor(eta, dv[, "bs"]),
      cor_logndt_bs = stats::cor(log(nd), dv[, "bs"]),
      stringsAsFactors = FALSE)
    cat(".")
    utils::flush.console()
  }
  cat(" nt=", nt, "\n", sep = "")
}
d <- do.call(rbind, rows)

cat("\n== the nuisance term against the trial count ==\n")
cat("paired on learners: the same seeds give the same hundred learners",
    " at every count\n", sep = "")
agg <- do.call(rbind, lapply(counts, function(nt) {
  s <- d[d$nt == nt, , drop = FALSE]
  data.frame(nt = nt, n = nrow(s),
             mean_log_m = mean(s$mean_log_m),
             sd_log_m = mean(s$sd_log_m),
             sd_log_ndt = mean(s$sd_log_ndt),
             sd_eta = mean(s$sd_eta),
             share_min_pct = 100 * mean(s$share_min),
             cor_eta_bs = mean(s$cor_eta_bs),
             cor_logndt_bs = mean(s$cor_logndt_bs),
             stringsAsFactors = FALSE)
}))
agg[-(1:2)] <- lapply(agg[-(1:2)], function(x) round(x, 4))
print(agg, row.names = FALSE)

cat("\n== paired differences, 400 trials minus 100 ==\n")
for (k in c("mean_log_m", "sd_log_m", "sd_eta", "share_min",
            "cor_eta_bs")) {
  x1 <- d[[k]][d$nt == 100L]
  x4 <- d[[k]][d$nt == 400L]
  tt <- stats::t.test(x4, x1, paired = TRUE)
  cat(sprintf("  %-12s %+8.4f  95%% CI (%+.4f, %+.4f)  p = %s\n", k,
              mean(x4 - x1), tt$conf.int[[1L]], tt$conf.int[[2L]],
              format(tt$p.value, digits = 3)))
}
cat("\nwith ", NS, " paired seeds a difference in sd(log m) smaller",
    " than about ", round(2 * stats::sd(d$sd_log_m[d$nt == 200L]) /
                            sqrt(NS), 4),
    " is not resolvable here.\n", sep = "")
saveRDS(d, "dev/learnhier-mvsn.rds")
cat("wrote dev/learnhier-mvsn.rds\n")
