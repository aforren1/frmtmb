# Inner-rule variants on the weighted Bernoulli and the mu-sigma case.
#   NQ_VARIANT=both|mode|scale Rscript dev/nestquad-variants.R
# linear = the package build; the variants come from the saved adapted
# build (dev/nestquad-adapt-version.R.txt). Seed 20261020.
.libPaths(c("C:/Users/adf44/source/r/wt-nestquad-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-nestquad"
source(file.path(wt, "dev/nestquad-ref.R"))
ad <- new.env(parent = ns)
sys.source(file.path(wt, "dev/nestquad-adapt-version.R.txt"), envir = ad)
set.seed(20261020)
ng <- 12; nt <- 6; n <- ng * nt
d <- data.frame(ga = factor(rep(seq_len(ng), each = nt)),
                gb = factor(rep(c("i", "ii"), length.out = n)),
                x = rnorm(n), w = runif(n, 0.5, 2))
ab <- interaction(d$ga, d$gb, drop = TRUE)
d$eta <- 0.2 + 0.5 * d$x + rnorm(ng, 0, 0.7)[d$ga] +
  rnorm(nlevels(ab), 0, 1)[ab]
d$y <- rbinom(n, 1, plogis(d$eta))
fit <- frm(bf(y | weights(w) ~ 1 + x + (1 | ga/gb)), family = bernoulli(),
           data = d, quadrature = TRUE)
fr <- fit$frame; par <- fit$opt$par
pl <- ns$imp_par_list(fr$par_template, par, fr$map)
th <- pl$theta
s <- vapply(fr$re_blocks, function(bk)
  sqrt(ns$cs_entry(bk$covstruct)$vcov(th[bk$theta_idx], bk)[1, 1]), 1)
b <- unname(pl$beta)
ref <- ref_nested_gh("bernoulli", d$y, b[1] + b[2] * d$x, d$ga, ab, s[1],
                     s[2], w = d$w, k = 120)
lap <- RTMB::MakeADFun(ns$build_objective(fr), fr$par_template, random = "b",
                       map = fr$map, silent = TRUE)
plan <- ns$quad_check_blocks(fr)
anc <- ns$quad_nest_anchor(lap, plan, par)
v <- Sys.getenv("NQ_VARIANT")
for (kk in list(c(11, 11), c(11, 15), c(15, 15), c(21, 21), c(31, 31))) {
  for (who in c("linear", v)) {
    obj_fun <- if (who == "linear") ns$quad_nest_objective else
      ad$quad_nest_objective
    qo <- obj_fun(fr, plan, anc, kk, pl)
    qo$prime(pl)
    val <- -qo$fn(pl)
    cat(sprintf("BERN %s ko %d ki %d err %.3e\n", who, kk[1], kk[2],
                abs(val - ref)))
  }
}
