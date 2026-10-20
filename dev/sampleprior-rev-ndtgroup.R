# Reviewer of lane sampleprior, fix round 1: wiener() with ndt_group(),
# where ndt is a FRACTION of each subject's fastest response; what does
# the default "uniform(0, min_Y)" (and its bound) do there?
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
suppressPackageStartupMessages(library(frmtmb.eam))
set.seed(20261120)
n <- 120
d <- frmtmb.eam::ddm_simulate(n, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
d$subj <- factor(rep(1:4, length.out = n))
d$x <- rnorm(n)
cat("min rt", min(d$rt), " per subject", tapply(d$rt, d$subj, min), "\n")
for (f in list(rt | dec(upper) ~ x, rt | dec(upper) + ndt_group(subj) ~ x)) {
  cat("\n==", deparse(f), "\n")
  r <- tryCatch({
    print(default_prior(f, data = d, family = frmtmb.eam::wiener(),
                        route = "sample"))
    od <- cap_obj(f, d, frmtmb.eam::wiener())
    cat("par:", names(od$par), "\nlower:", CAP$args$lower,
        "\nupper:", CAP$args$upper, "\n")
    fit <- suppressMessages(suppressWarnings(asNamespace("frmtmb.sample")$
      sample_assemble(f, d, frmtmb.eam::wiener(), data2 = list(),
                      start = NULL, fit_control = frmtmb_control(),
                      na.action = stats::na.omit, REML = FALSE)))
    lp <- Filter(function(l) identical(l$dpar, "ndt"), fit$frame$linpreds)[[1]]
    cat("ndt link", lp$link$name, " image at -30, 30:",
        lp$link$linkinv(c(-30, 30)), "\n")
    "ok"
  }, error = function(e) paste("ERROR", conditionMessage(e)))
  cat(r, "\n")
}
