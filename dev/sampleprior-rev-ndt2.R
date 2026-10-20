# Reviewer of lane sampleprior, fix round 2: the ndt default and its bound
# for wiener(), lba() and rdm(), with and without ndt_group(); the
# ungrouped wiener against brms 2.23.0's default_prior().
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
suppressPackageStartupMessages(library(frmtmb.eam))
set.seed(20261120)
dw <- frmtmb.eam::ddm_simulate(120, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
dw$subj <- factor(rep(1:4, length.out = 120)); dw$x <- rnorm(120)
dl <- frmtmb.eam::lba_simulate(200, v = c(2.4, 1.4), A = 0.5, k = 0.4, ndt = 0.2)
dl$subj <- factor(rep(1:4, length.out = 200))
dr <- frmtmb.eam::rdm_simulate(200, v = c(2.4, 1.4), A = 0.5, k = 0.5, ndt = 0.2)
dr$subj <- factor(rep(1:4, length.out = 200))
cases <- list(
  list("wiener", rt | dec(upper) ~ x, dw, frmtmb.eam::wiener()),
  list("wiener + ndt_group", rt | dec(upper) + ndt_group(subj) ~ x, dw,
       frmtmb.eam::wiener()),
  list("lba", rt | vint(choice) ~ 1, dl, frmtmb.eam::lba(2)),
  list("lba + ndt_group", rt | vint(choice) + ndt_group(subj) ~ 1, dl,
       frmtmb.eam::lba(2)),
  list("rdm", rt | vint(choice) ~ 1, dr, frmtmb.eam::rdm(2)),
  list("rdm + ndt_group", rt | vint(choice) + ndt_group(subj) ~ 1, dr,
       frmtmb.eam::rdm(2)))
for (cs in cases) {
  r <- tryCatch({
    tab <- as.data.frame(default_prior(cs[[2]], data = cs[[3]],
                                       family = cs[[4]], route = "sample"))
    nd <- tab$prior[tab$class == "ndt"]
    od <- cap_obj(cs[[2]], cs[[3]], cs[[4]])
    up <- CAP$args$upper
    k <- if (is.null(up)) integer(0) else which(is.finite(up))
    sprintf("ndt default %s | min rt %.15g | finite upper bounds: %s",
            paste(nd, collapse = ","), min(cs[[3]]$rt),
            if (length(k)) paste(names(od$par)[k], up[k], collapse = " ") else "none")
  }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 100)))
  cat(sprintf("%-20s %s\n", cs[[1]], r))
}
bp <- as.data.frame(brms::default_prior(rt | dec(upper) ~ x, data = dw,
                                        family = brms::wiener()))
cat("brms wiener ndt row:", bp$prior[bp$class == "ndt"], " ub", bp$ub[bp$class == "ndt"],
    " (min_Y =", sprintf("%.15g", min(dw$rt)), ")\n")
# the posterior under ndt_group(): must cover the ML fraction
fit <- suppressWarnings(frm(rt | dec(upper) + ndt_group(subj) ~ x, data = dw,
                            family = frmtmb.eam::wiener()))
ml <- frm_linpred(fit, dpar = "ndt", type = "response")[1]
ds <- suppressMessages(suppressWarnings(frm_sample(
  fit, chains = 2, iter = 600, refresh = 0, seed = 4)))
v <- ds$draws[, "ndt"]
cat(sprintf("wiener + ndt_group, defaults: ML fraction %.4f | posterior q01/05/50/95/99 %s | rhat %.3f\n",
            ml, paste(sprintf("%.4f", quantile(v, c(.01, .05, .5, .95, .99))), collapse = " "),
            rstan::Rhat(matrix(v, ncol = 2))))
