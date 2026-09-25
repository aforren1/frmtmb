# Items 3.4 and 3.5: recovery and Wald coverage at 30 subjects by 400
# trials, the realistic eam scale in dev/extension-gaps-plan.md.
#
# Usage: Rscript dev/phase3b-eam-recovery.R <arm> <seed_from> <seed_to>
#   arm = cens      deadline censoring at 1.0 s, no contamination
#         cont      5 percent uniform contaminants over [0.1, 5] s,
#                   contaminant = TRUE and max_ndt = 0.5
#         contdef   as cont, fitted at the default ndt bound
#         collapse  NO contaminants, fitted with contaminant = TRUE
# One RDS per seed in dev/phase3b-log/recov/<arm>-<seed>.rds, written
# only when the fit returns, so counts are taken from the files.
#
# Truth, shared by all arms: mu = 0.4 + 0.9 * cond + u_s, u_s ~ N(0,
# 0.35); log bs = log 1.4 + b_s, b_s ~ N(0, 0.20); ndt = 0.25; bias 0.5.
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
args <- commandArgs(trailingOnly = TRUE)
arm <- args[[1]]
seeds <- seq(as.integer(args[[2]]), as.integer(args[[3]]))
stopifnot(arm %in% c("cens", "cont", "contdef", "collapse"))
dir.create("dev/phase3b-log/recov", showWarnings = FALSE, recursive = TRUE)

NS <- 30L; NT <- 400L
DEADLINE <- 1.0
LAMBDA <- 0.05
CLO <- 0.1; CHI <- 5
MAXNDT <- 0.5
truth <- c(mu_int = 0.4, mu_cond = 0.9, bs_int = log(1.4),
           ndt = 0.25, sd_mu = 0.35, sd_bs = 0.20,
           lambda = qlogis(LAMBDA))

simulate_one <- function(seed) {
  set.seed(seed)
  u <- rnorm(NS, 0, 0.35)
  b <- rnorm(NS, 0, 0.20)
  out <- vector("list", NS)
  for (s in seq_len(NS)) {
    cond <- rep(0:1, length.out = NT)
    mu <- 0.4 + 0.9 * cond + u[s]
    dd <- ddm_simulate(NT, mu = mu, bs = 1.4 * exp(b[s]), ndt = 0.25,
                       bias = 0.5)
    dd$cond <- cond
    dd$s <- factor(s, levels = seq_len(NS))
    out[[s]] <- dd
  }
  d <- do.call(rbind, out)
  if (arm %in% c("cont", "contdef")) {
    hit <- runif(nrow(d)) < LAMBDA
    d$rt[hit] <- runif(sum(hit), CLO, CHI)
    d$upper[hit] <- rbinom(sum(hit), 1, 0.5)
    attr(d, "n_cont") <- sum(hit)
  }
  if (arm == "cens") {
    late <- d$rt > DEADLINE
    d$rt[late] <- DEADLINE
    d$cens <- ifelse(late, "right", "none")
    # the boundary of a censored trial is not observed; the family must
    # not read it, so it is set to a value that would be wrong half the
    # time if it were read
    d$upper[late] <- 0
  }
  d
}

fit_one <- function(d) {
  if (arm == "cens") {
    f <- bf(rt | dec(upper) + cens(cens) ~ cond + (1 | s),
            bs ~ 1 + (1 | s), ndt ~ 1, bias = 0.5)
    fam <- wiener()
  } else {
    f <- bf(rt | dec(upper) ~ cond + (1 | s), bs ~ 1 + (1 | s), ndt ~ 1,
            lambda ~ 1, bias = 0.5)
    # `cont` states a bound above the fastest response, which the
    # contaminant makes a model rather than a mistake: some guesses are
    # faster than the true non-decision time. `contdef` keeps the
    # default bound on the same draws, to measure what that costs.
    fam <- if (arm == "cont") {
      wiener(contaminant = TRUE, max_ndt = MAXNDT)
    } else {
      wiener(contaminant = TRUE)
    }
  }
  frm(f, family = fam, data = d)
}

for (seed in seeds) {
  t0 <- proc.time()[["elapsed"]]
  d <- simulate_one(seed)
  res <- tryCatch({
    fit <- fit_one(d)
    fe <- unlist(fixef_by_dpar(fit))
    se <- tryCatch(sqrt(diag(as.matrix(vcov(fit)))),
                   error = function(e) NULL)
    ci <- tryCatch(confint(fit), error = function(e) NULL)
    diag_out <- tryCatch(utils::capture.output(diagnose(fit)),
                         error = function(e) conditionMessage(e))
    list(seed = seed, arm = arm, fixef = fe, se = se, confint = ci,
         varcorr = tryCatch(VarCorr(fit), error = function(e) NULL),
         loglik = as.numeric(logLik(fit)),
         ndt_ub = frmtmb::single_response(fit, "fit")$family$ndt_bound$ub,
         crange = frmtmb::single_response(fit, "fit")$family$contaminant_range,
         conv = fit$opt$convergence,
         pdhess = tryCatch(fit$sdr$pdHess, error = function(e) NA),
         diagnose = diag_out,
         n_cens = if (arm == "cens") sum(d$cens == "right") else NA,
         n_cont = if (arm %in% c("cont", "contdef")) attr(d, "n_cont") else 0L,
         elapsed = proc.time()[["elapsed"]] - t0,
         pkg = c(frmtmb = as.character(packageVersion("frmtmb")),
                 eam = as.character(packageVersion("frmtmb.eam"))))
  }, error = function(e) {
    list(seed = seed, arm = arm, error = conditionMessage(e))
  })
  if (arm == "collapse" && is.null(res$error)) {
    # the plain model on the same data, for the likelihood comparison
    fp <- tryCatch(frm(bf(rt | dec(upper) ~ cond + (1 | s),
                          bs ~ 1 + (1 | s), ndt ~ 1, bias = 0.5),
                       family = wiener(), data = d),
                   error = function(e) NULL)
    res$loglik_plain <- if (is.null(fp)) NA else as.numeric(logLik(fp))
  }
  saveRDS(res, sprintf("dev/phase3b-log/recov/%s-%d.rds", arm, seed))
  cat(arm, seed, if (is.null(res$error)) "ok" else res$error,
      round(proc.time()[["elapsed"]] - t0, 1), "\n")
}
