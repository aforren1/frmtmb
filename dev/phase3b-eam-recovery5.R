# Punch round 1: recovery and Wald coverage at 30 subjects by 400
# trials on the build with B1, M1 and decisions (a) and (b).
#
# Usage: Rscript dev/phase3b-eam-recovery5.R <arm> <seed_from> <seed_to>
#   cont      5 percent contaminants uniform on [0.1, 5] s, coin-flip
#             boundary; fitted with contaminant_range = range(rt), the
#             old default, now only reachable by passing it; max_ndt 0.5
#   contfix   the same draws, contaminant_range = c(0.1, 5)
#   contdl    the same draws with a 3 s deadline: rows slower than 3 s
#             are not recorded, trunc(ub = 3) declares it, and the
#             window is the default under a deadline
#   collapse  NO contaminants, contaminant_range = c(0, 5)
#   cleft     no contaminants; rows faster than 0.45 s are recorded as
#             LEFT-censored at 0.45 s, and every fifth other row as
#             INTERVAL-censored into its 100 ms bin; dec() is the true
#             boundary on both, and is read
#   cens      no contaminants; rows slower than 1.0 s right-censored at
#             1.0 s with dec() set to 0, unread (the first round's arm,
#             rerun on the build with the rewritten survival)
# Writes dev/phase3b-log/recov3/<arm>-<seed>.rds, the fit or its error.
#
# Truth: mu = 0.4 + 0.9 cond + u_s, u_s ~ N(0, 0.35); log bs = log 1.4 +
# b_s, b_s ~ N(0, 0.20); ndt = 0.25; bias 0.5.
.libPaths(c(Sys.getenv("P3B_LIB", "C:/Users/adf44/source/r/phase3b-lib"),
            "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
args <- commandArgs(trailingOnly = TRUE)
arm <- args[[1]]
seeds <- seq(as.integer(args[[2]]), as.integer(args[[3]]))
stopifnot(arm %in% c("cont", "contfix", "contdl", "collapse", "cleft", "cens"))
dir.create("dev/phase3b-log/recov3", showWarnings = FALSE, recursive = TRUE)
NS <- 30L; NT <- 400L

simulate_one <- function(seed) {
  set.seed(seed)
  u <- rnorm(NS, 0, 0.35)
  b <- rnorm(NS, 0, 0.20)
  d <- do.call(rbind, lapply(seq_len(NS), function(s) {
    cond <- rep(0:1, length.out = NT)
    x <- ddm_simulate(NT, mu = 0.4 + 0.9 * cond + u[s],
                      bs = 1.4 * exp(b[s]), ndt = 0.25, bias = 0.5)
    x$cond <- cond
    x$s <- factor(s, levels = seq_len(NS))
    x
  }))
  if (arm %in% c("cont", "contfix", "contdl")) {
    hit <- runif(nrow(d)) < 0.05
    d$rt[hit] <- runif(sum(hit), 0.1, 5)
    d$upper[hit] <- rbinom(sum(hit), 1, 0.5)
    d$cont <- hit
  }
  if (arm == "contdl") d <- d[d$rt < 3, ]
  if (arm == "cens") {
    late <- d$rt > 1.0
    d$rt[late] <- 1.0
    d$code <- as.integer(late)
    d$upper[late] <- 0
  }
  if (arm == "cleft") {
    d$code <- 0L
    d$y2 <- d$rt
    fast <- d$rt < 0.45
    d$code[fast] <- -1L
    d$rt[fast] <- 0.45
    pick <- which(!fast)[seq(5, sum(!fast), by = 5)]
    lo <- floor(d$rt[pick] * 10) / 10
    d$code[pick] <- 2L
    d$y2[pick] <- lo + 0.1
    d$rt[pick] <- pmax(lo, 0.45)
  }
  d
}

fit_one <- function(d) {
  base <- rt | dec(upper) ~ cond + (1 | s)
  if (arm == "cens") {
    return(frm(bf(rt | dec(upper) + cens(code) ~ cond + (1 | s),
                  bs ~ 1 + (1 | s), ndt ~ 1, bias = 0.5),
               family = wiener(), data = d))
  }
  if (arm == "cleft") {
    return(frm(bf(rt | dec(upper) + cens(code, y2) ~ cond + (1 | s),
                  bs ~ 1 + (1 | s), ndt ~ 1, bias = 0.5),
               family = wiener(), data = d))
  }
  if (arm == "collapse") {
    return(frm(bf(base, bs ~ 1 + (1 | s), ndt ~ 1, lambda ~ 1, bias = 0.5),
               family = wiener(contaminant = TRUE,
                               contaminant_range = c(0, 5)), data = d))
  }
  if (arm == "contdl") {
    return(frm(bf(rt | dec(upper) + trunc(ub = 3) ~ cond + (1 | s),
                  bs ~ 1 + (1 | s), ndt ~ 1, lambda ~ 1, bias = 0.5),
               family = wiener(contaminant = TRUE, max_ndt = 0.5),
               data = d))
  }
  rng <- if (arm == "contfix") c(0.1, 5) else range(d$rt)
  frm(bf(base, bs ~ 1 + (1 | s), ndt ~ 1, lambda ~ 1, bias = 0.5),
      family = wiener(contaminant = TRUE, max_ndt = 0.5,
                      contaminant_range = rng), data = d)
}

for (seed in seeds) {
  t0 <- proc.time()[["elapsed"]]
  d <- simulate_one(seed)
  res <- tryCatch({
    fit <- fit_one(d)
    fam <- stats::family(fit)
    list(seed = seed, arm = arm, fixef = unlist(fixef_by_dpar(fit)),
         confint = tryCatch(confint(fit), error = function(e)
           conditionMessage(e)),
         loglik = as.numeric(logLik(fit)),
         ndt_ub = fam$ndt_bound$ub, crange = fam$contaminant_range,
         conv = fit$opt$convergence,
         diagnose = tryCatch(utils::capture.output(diagnose(fit)),
                             error = function(e) conditionMessage(e)),
         n = nrow(d),
         n_cont = if (is.null(d$cont)) 0L else sum(d$cont),
         n_left = if (is.null(d$code)) 0L else sum(d$code == -1L),
         n_int = if (is.null(d$code)) 0L else sum(d$code == 2L),
         elapsed = proc.time()[["elapsed"]] - t0,
         eam_path = find.package("frmtmb.eam"))
  }, error = function(e) list(seed = seed, arm = arm,
                              error = conditionMessage(e),
                              eam_path = find.package("frmtmb.eam")))
  if (arm == "collapse" && is.null(res$error)) {
    fp <- tryCatch(frm(bf(rt | dec(upper) ~ cond + (1 | s),
                          bs ~ 1 + (1 | s), ndt ~ 1, bias = 0.5),
                       family = wiener(), data = d),
                   error = function(e) NULL)
    res$loglik_plain <- if (is.null(fp)) NA else as.numeric(logLik(fp))
  }
  saveRDS(res, sprintf("dev/phase3b-log/recov3/%s-%d.rds", arm, seed))
  cat(arm, seed, if (is.null(res$error)) "ok" else res$error,
      round(proc.time()[["elapsed"]] - t0, 1), "\n")
}
