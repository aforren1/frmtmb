# Rerun of dev/phase3b-cont-window.R on the final build, into contwin2/.
# User decision (a), 2026-09-24: which window should a contaminant
# spread over when the user gives none? Measured on one subject of 4000
# trials, 25 seeds, mu 0.5, bs 2.0, ndt 0.3, bias 0.5, contaminant share
# 0.05 uniform on [0, 5] s with a coin-flip boundary.
#   deadline: rows slower than 5 s are not recorded, and the model
#     declares it with trunc(ub = 5). Windows: observed [min, max],
#     [min, 5] (the trunc bound as the top), and [0, 5].
#   none: nothing is dropped and no trunc(). Windows: observed
#     [min, max] and the true [0, 5].
# The truth for lambda is the share among RECORDED rows, which a
# deadline raises slightly by dropping slow diffusion trials.
# Usage: Rscript dev/phase3b-cont-window.R <design> <seed_from> <seed_to>
# Output: dev/phase3b-log/contwin2/<design>-<seed>.rds
.libPaths(c(Sys.getenv("P3B_LIB", "C:/Users/adf44/source/r/phase3b-lib"),
            "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
a <- commandArgs(trailingOnly = TRUE)
design <- a[[1]]
dir.create("dev/phase3b-log/contwin2", showWarnings = FALSE)
for (seed in seq(as.integer(a[[2]]), as.integer(a[[3]]))) {
  set.seed(seed)
  n <- 4000
  d <- ddm_simulate(n, mu = 0.5, bs = 2.0, ndt = 0.3, bias = 0.5)
  hit <- runif(n) < 0.05
  d$rt[hit] <- runif(sum(hit), 0, 5)
  d$upper[hit] <- rbinom(sum(hit), 1, 0.5)
  d$cont <- hit
  if (design == "deadline") d <- d[d$rt < 5, ]
  truth <- mean(d$cont)
  wins <- if (design == "deadline") {
    list(observed = NULL, top5 = c(min(d$rt), 5), zero5 = c(0, 5))
  } else {
    list(observed = range(d$rt), zero5 = c(0, 5))
  }
  f <- if (design == "deadline") {
    bf(rt | dec(upper) + trunc(ub = 5) ~ 1, bias = 0.5)
  } else {
    bf(rt | dec(upper) ~ 1, bias = 0.5)
  }
  res <- lapply(names(wins), function(nm) {
    fit <- tryCatch(suppressWarnings(frm(
      f, family = wiener(contaminant = TRUE, max_ndt = 0.6,
                         contaminant_range = wins[[nm]]), data = d)),
      error = function(e) e)
    if (inherits(fit, "error")) return(list(window = nm,
                                            error = conditionMessage(fit)))
    ci <- confint(fit)["lambda_(Intercept)", ]
    list(window = nm, est = ci[["est"]], lwr = ci[["lwr"]], upr = ci[["upr"]],
         crange = stats::family(fit)$contaminant_range,
         conv = fit$opt$convergence)
  })
  saveRDS(list(seed = seed, design = design, truth = truth, n = nrow(d),
               res = res, eam = find.package("frmtmb.eam")),
          sprintf("dev/phase3b-log/contwin2/%s-%d.rds", design, seed))
  cat(design, seed, "\n")
}
