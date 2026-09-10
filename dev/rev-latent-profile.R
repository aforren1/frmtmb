# Reviewer, punch round 2, item 1 second half: does
# `confint(method = "profile")` actually recover nominal coverage on the
# two coefficients `?lca` now sends users there for?
#
# THE OBSTACLE, AND HOW IT IS HANDLED. The quantity the recovery table
# scores is a CONTRAST, `theta_{perm[a]} - theta_{perm[1]}`, because
# frmtmb references the last class and the truth references class 1.
# `confint(method = "profile")` profiles a single named parameter, not a
# contrast, so it applies to the scored quantity only where
# `perm[1] == K`: there `hit(K, -1)` contributes nothing and the
# contrast IS the raw parameter `theta_{perm[a]}_x2`.
#
# On the fresh block that is 59 of the 200 replicates. Those 59 are the
# subset used here, and the Wald coverage is recomputed on the SAME
# subset so the comparison is like for like: a profile interval that
# beat Wald on a different subset would prove nothing.
#
# Seeds: the fresh block's own, 20270401 to 20270600, restricted to the
# 59 whose `perm_truth` starts with 4. Script writes
# dev/rev-latent-profile.tsv.
#
#   Rscript dev/rev-latent-profile.R [nrep]

source("C:/Users/adf44/source/r/frmtmb-wt-latent/dev/rev-latent-env3.R")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.latent)
})
D <- "C:/Users/adf44/source/r/frmtmb-wt-latent/dev/"
source(paste0(D, "latent-lca-sim.R"))
K <- 4L
tight <- frmtmb_control(optCtrl = list(iter.max = 20000, eval.max = 20000,
                                       rel.tol = 1e-14, x.tol = 1e-14))
truth_gam <- lca_truth_gam          # rows: class 1..4, cols: int, x1, x2

r <- utils::read.delim(paste0(D, "latent-2p4-recheck-oos.tsv"),
                       stringsAsFactors = FALSE)
r$perm_truth <- trimws(as.character(r$perm_truth))
sub <- r[substr(r$perm_truth, 1L, 1L) == "4", ]
args <- commandArgs(trailingOnly = TRUE)
NR <- if (length(args)) min(as.integer(args[[1L]]), nrow(sub)) else nrow(sub)
sub <- sub[seq_len(NR), ]
cat("fresh-block replicates whose perm starts with 4:", nrow(sub), "\n")
cat("on those the scored contrast IS a raw parameter.\n\n")

out <- paste0(D, "rev-latent-profile.tsv")
cat(paste(c("seed", "perm", "coef", "par", "est", "truth",
            "wald_lo", "wald_hi", "prof_lo", "prof_hi",
            "wald_cov", "prof_cov", "wald_w", "prof_w", "secs"),
          collapse = "\t"), "\n", sep = "", file = out)

q <- stats::qnorm(0.975)
t0all <- Sys.time()
for (i in seq_len(nrow(sub))) {
  seed <- sub$seed[i]
  perm <- as.integer(strsplit(sub$perm_truth[i], "")[[1L]])
  s <- lca_sim(seed = seed)
  fit <- suppressWarnings(suppressMessages(
    frm(bf(Y ~ x1 + x2), family = lca(K = K), data = s$dd,
        control = tight)))
  V <- vcov(fit)
  b <- fixef(fit)
  cn <- names(b[[1L]])[3L]                   # the binary covariate x2
  for (a in c(3L, 4L)) {                     # truth classes 3 and 4
    k <- perm[a]
    stopifnot(k < K)                         # else it is not a raw par
    pnm <- paste0("theta", k, "_", cn)
    est <- unname(b[[paste0("theta", k)]][3L])
    se <- sqrt(V[pnm, pnm])
    tv <- truth_gam[a, 3L]
    t1 <- Sys.time()
    pc <- try(suppressWarnings(stats::confint(fit, parm = pnm,
                                              method = "profile")),
              silent = TRUE)
    secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    plo <- if (inherits(pc, "try-error")) NA_real_ else pc[1L, "lwr"]
    phi <- if (inherits(pc, "try-error")) NA_real_ else pc[1L, "upr"]
    wlo <- est - q * se; whi <- est + q * se
    cat(paste(c(seed, sub$perm_truth[i], paste0("c", a, "_3"), pnm,
                formatC(c(est, tv, wlo, whi, plo, phi), digits = 6,
                        format = "g"),
                as.integer(tv > wlo & tv < whi),
                if (is.na(plo)) NA_integer_ else
                  as.integer(tv > plo & tv < phi),
                formatC(c(whi - wlo, phi - plo, secs), digits = 5,
                        format = "g")),
              collapse = "\t"), "\n", sep = "", file = out, append = TRUE)
  }
  if (i %% 10L == 0L) {
    cat("  ", i, " of ", nrow(sub), "  (",
        format(as.numeric(difftime(Sys.time(), t0all, units = "mins")),
               digits = 3), " min)\n", sep = "")
    flush(stdout())
  }
}

## ------------------------------------------------------------ summary
d <- utils::read.delim(out)
mci <- function(k, n) {
  p <- k / n
  se <- sqrt(p * (1 - p) / n)
  sprintf("%.4f (%d/%d)  MC 95%%: %.4f to %.4f", p, k, n,
          p - 1.96 * se, p + 1.96 * se)
}
cat("\n==== the 59-replicate subset, Wald against profile ====\n")
cat("profile intervals that failed to compute:", sum(is.na(d$prof_cov)),
    "of", nrow(d), "\n")
ok <- !is.na(d$prof_cov)
for (cf in c("c3_3", "c4_3")) {
  z <- d[ok & d$coef == cf, ]
  cat(sprintf("\n%s (%s)\n", cf, unique(z$par)[1L]))
  cat("  Wald    :", mci(sum(z$wald_cov), nrow(z)), "\n")
  cat("  profile :", mci(sum(z$prof_cov), nrow(z)), "\n")
  cat("  median interval width, Wald", format(stats::median(z$wald_w),
                                              digits = 4),
      " profile", format(stats::median(z$prof_w), digits = 4),
      " ratio", format(stats::median(z$prof_w / z$wald_w), digits = 4),
      "\n")
}
z <- d[ok, ]
cat("\nboth coefficients pooled on this subset:\n")
cat("  Wald    :", mci(sum(z$wald_cov), nrow(z)), "\n")
cat("  profile :", mci(sum(z$prof_cov), nrow(z)), "\n")
cat("  profile wider than Wald on", sum(z$prof_w > z$wald_w), "of",
    nrow(z), "intervals\n")
cat("  median width ratio, profile over Wald:",
    format(stats::median(z$prof_w / z$wald_w), digits = 4), "\n")
cat("  median seconds per profile interval:",
    format(stats::median(z$secs), digits = 4), "\n")
cat("  total wall clock:",
    format(as.numeric(difftime(Sys.time(), t0all, units = "mins")),
           digits = 4), "minutes\n")
