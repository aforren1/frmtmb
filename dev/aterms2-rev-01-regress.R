# Reviewer, claim 1 (no regression). Models WITHOUT the new terms that
# the per-response frame, the Z row count, the cs() matrix, the mi()
# loop and the rewritten count densities pass through. Each output is
# stored whole and compared with identical() by
# dev/aterms2-rev-01-compare.R. Run once per arm, one process each:
#   Rscript dev/aterms2-rev-01-regress.R base
#   Rscript dev/aterms2-rev-01-regress.R lane
# Seeds: data 101, perturbation 7, predict 11, simulate 13.
arm <- commandArgs(TRUE)[[1]]
libs <- c("C:/Users/adf44/source/r/wt-aterms2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

set.seed(101)
n <- 120
G <- 8
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                g = factor(rep(seq_len(G), length.out = n)),
                f = factor(rep(c("a", "b"), each = n / 2)),
                gg = factor(rep(c("u", "v"), n / 2)),
                t = rep(seq_len(n / G), each = G),
                wt = runif(n, 0.5, 2), nt = sample(5:15, n, TRUE),
                time = runif(n, 0.5, 3))
re <- rnorm(G, 0, 0.6)
d$fg <- factor(ifelse(as.integer(d$g) <= G / 2, "a", "b"))
d$y1 <- 1 + d$x + re[d$g] + rnorm(n)
d$y2 <- 0.5 - d$z + 0.5 * re[d$g] + rnorm(n)
d$yc <- rpois(n, exp(0.3 + 0.4 * d$x + re[d$g]))
d$yn <- rnbinom(n, mu = exp(0.5 + 0.4 * d$x), size = 2)
d$yz <- ifelse(runif(n) < 0.3, 0L, d$yc)
d$yb <- rbinom(n, d$nt, plogis(-0.2 + 0.5 * d$x))
cut4 <- function(v) factor(cut(v, c(-Inf, -0.5, 0.3, 1.2, Inf),
                               labels = FALSE), ordered = TRUE)
d$o1 <- cut4(d$x + rnorm(n))
d$o2 <- cut4(d$z + rnorm(n))
d$om <- factor(sample(1:4, n, TRUE), ordered = TRUE)
d$xm <- d$w + rnorm(n, 0, 0.5)
d$ym <- 1 + 0.7 * d$xm + 0.3 * d$z + rnorm(n, 0, 0.5)
d$xm[c(3, 9, 40, 77)] <- NA
d$sdx <- runif(n, 0.1, 0.4)
d$xo <- d$x + rnorm(n, 0, d$sdx)
d$cc <- rep(c(0, 0, 1, -1), n / 4)
d$cc3 <- rep(c(0, 2, 0, 1), n / 4)
d$yup <- d$y1 + runif(n, 0.2, 1)
d$ytr <- pmax(d$y1, -1)
d$ypt <- pmax(d$yc, 1L)
d$y2na <- d$y2
d$y2na[c(5, 17, 60)] <- NA
d$y1na <- d$y1
d$y1na[c(2, 50)] <- NA
d$xna <- d$x
d$xna[c(8, 90)] <- NA

nd <- d[c(1:6, 60:64), ]
fams <- list(
  mv_rescor_gauss = list(
    function() frm(bf(y1 ~ x) + bf(y2 ~ z) + set_rescor(TRUE), data = d,
                   family = gaussian())),
  mv_rescor_student = list(
    function() frm(bf(y1 ~ x) + bf(y2 ~ z) + set_rescor(TRUE), data = d,
                   family = student())),
  mv_ordinal = list(
    function() frm(bf(o1 ~ x) + bf(o2 ~ z + (1 | g)), data = d,
                   family = cumulative())),
  mv_ordinal_cs = list(
    function() frm(bf(o1 ~ cs(x)) + bf(o2 ~ z), data = d,
                   family = sratio())),
  mv_ordinal_thronly = list(
    function() frm(bf(o1 ~ 1) + bf(o2 ~ 1 + (1 | g)), data = d,
                   family = cumulative())),
  mv_id_mixed = list(
    function() frm(bf(y1 ~ x + (1 + x | p | g)) + gaussian() +
                     bf(yc ~ z + (1 | p | g)) + poisson(), data = d)),
  mv_mi = list(
    function() frm(bf(ym ~ mi(xm) + z) + bf(xm | mi() ~ w), data = d,
                   family = gaussian())),
  mv_mi_inter = list(
    function() frm(bf(ym ~ mi(xm) * z) + bf(xm | mi() ~ w), data = d,
                   family = gaussian())),
  me = list(function() frm(y1 ~ me(xo, sdx) + z, data = d)),
  cens_gauss = list(function() frm(y1 | cens(cc) ~ x, data = d)),
  cens_interval = list(function() frm(y1 | cens(cc3, yup) ~ x, data = d)),
  cens_pois = list(function() frm(yc | cens(cc) ~ x, data = d,
                                  family = poisson())),
  trunc_gauss = list(function() frm(ytr | trunc(lb = -1) ~ x, data = d)),
  trunc_pois = list(function() frm(ypt | trunc(lb = 1) ~ x, data = d,
                                   family = poisson())),
  weights_re = list(function() frm(y1 | weights(wt) ~ x + (1 | g),
                                   data = d)),
  trials = list(function() frm(yb | trials(nt) ~ x + (1 | g), data = d,
                               family = binomial())),
  thres_gr = list(function() frm(o1 | thres(gr = gg) ~ x, data = d,
                                 family = cumulative())),
  gr_by = list(function() frm(y1 ~ x + (1 | gr(g, by = fg)), data = d)),
  mv_ar_percresp = list(
    function() frm(bf(y1 ~ x + ar(t, g)) + bf(y2 ~ z), data = d,
                   family = gaussian())),
  mv_na_one_resp = list(
    function() frm(bf(y1 ~ x) + bf(y2na ~ z), data = d,
                   family = gaussian())),
  mv_na_both_pred = list(
    function() frm(bf(y1na ~ xna) + bf(y2na ~ z + (1 | g)), data = d,
                   family = gaussian())),
  uni_na_omit = list(function() frm(y1na ~ xna + (1 | g), data = d)),
  uni_na_exclude = list(function() frm(y1na ~ xna + (1 | g), data = d,
                                       na.action = stats::na.exclude)),
  mv_na_exclude = list(
    function() frm(bf(y1na ~ x) + bf(y2 ~ z), data = d,
                   family = gaussian(), na.action = stats::na.exclude)),
  mv_smooth = list(
    function() frm(bf(y1 ~ s(x)) + bf(y2 ~ s(z)), data = d,
                   family = gaussian())),
  mv_mo = list(function() frm(bf(y1 ~ mo(om)) + bf(y2 ~ z), data = d,
                              family = gaussian())),
  pois_re = list(function() frm(yc ~ x + (1 | g), data = d,
                                family = poisson())),
  pois_identity = list(function() frm(yc ~ 1 + time, data = d,
                                      family = poisson("identity"))),
  pois_offset = list(function() frm(yc ~ x + offset(log(time)), data = d,
                                    family = poisson())),
  negbin = list(function() frm(bf(yn ~ x, shape ~ z), data = d,
                               family = negbinomial())),
  negbin_identity = list(function() frm(yn ~ 1 + time, data = d,
                                        family = negbinomial("identity"))),
  geometric = list(function() frm(yn ~ x, data = d, family = geometric())),
  geometric_sqrt = list(function() frm(yn ~ time, data = d,
                                       family = geometric("sqrt"))),
  zip = list(function() frm(yz ~ x, data = d,
                            family = zero_inflated_poisson())),
  zinb = list(function() frm(yz ~ x, data = d,
                             family = zero_inflated_negbinomial())),
  hurdle_nb = list(function() frm(yz ~ x, data = d,
                                  family = hurdle_negbinomial())),
  hurdle_pois = list(function() frm(yz ~ x, data = d,
                                    family = hurdle_poisson())),
  mix_pois = list(function() frm(yn ~ 1, data = d,
                                 family = mixture(poisson(), poisson()))),
  mv_poisson_mixed_links = list(
    function() frm(bf(yc ~ x) + poisson() + bf(yn ~ z) + negbinomial(),
                   data = d))
)

strip <- function(x) {
  if (is.function(x) || is.environment(x) || is.language(x)) return(NULL)
  if (methods::is(x, "Matrix")) return(as.matrix(x))
  if (is.data.frame(x)) x <- as.list(x)
  if (is.list(x)) {
    out <- lapply(x, strip)
    out <- out[!vapply(out, is.null, TRUE)]
    return(out)
  }
  if (is.atomic(x)) {
    at <- attributes(x)
    keep <- at[intersect(names(at), c("names", "dim", "dimnames",
                                      "levels", "class"))]
    attributes(x) <- keep
    return(x)
  }
  NULL
}
try_ <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)),
           error = function(e) paste("ERROR:", conditionMessage(e)))
}

res <- list()
for (nm in names(fams)) {
  cat(nm, "\n")
  f <- try_(fams[[nm]][[1]]())
  if (is.character(f)) {
    res[[nm]] <- list(fit = f)
    cat("  ", f, "\n")
    next
  }
  set.seed(7)
  p <- f$opt$par + rnorm(length(f$opt$par), 0, 0.05)
  r <- list()
  r$fn <- f$obj$fn(p)
  r$gr <- f$obj$gr(p)
  r$par <- f$opt$par
  r$logLik <- try_(logLik(f))
  r$nobs <- try_(nobs(f))
  r$fixef <- try_(fixef(f))
  r$coef <- try_(strip(coef(f)))
  r$ranef <- try_(strip(ranef(f)))
  r$VarCorr <- try_(strip(VarCorr(f)))
  r$frame <- strip(f$frame[setdiff(names(f$frame), "subset_rows")])
  resps <- names(f$spec$responses)
  for (rr in resps) {
    r[[paste0("fitted_", rr)]] <- try_(fitted(f, resp = rr))
    r[[paste0("fitted_nd_", rr)]] <- try_(fitted(f, newdata = nd, resp = rr))
    r[[paste0("linpred_", rr)]] <- try_(frm_linpred(f, resp = rr))
    set.seed(11)
    r[[paste0("predict_nd_", rr)]] <- try_(predict(f, newdata = nd,
                                                   resp = rr, ndraws = 20))
    r[[paste0("resid_", rr)]] <- try_(residuals(f, resp = rr))
  }
  r$fitted_all <- try_(fitted(f))
  r$sim <- try_(strip(simulate(f, nsim = 2, seed = 13)))
  set.seed(11)
  r$predict_all <- try_(predict(f, ndraws = 10))
  r$resid_pearson <- try_(residuals(f, type = "pearson"))
  r$ce <- if (nm %in% c("pois_re", "negbin", "mv_id_mixed", "mv_mi",
                        "trials", "cens_pois", "trunc_pois")) {
    try_(strip(lapply(conditional_effects(f), as.list)))
  }
  r$default_prior <- try_(as.data.frame(default_prior(f)))
  res[[nm]] <- r
}
out <- sprintf("C:/Users/adf44/source/r/frmtmb-wt-aterms2/dev/aterms2-rev-01-%s.rds",
               arm)
saveRDS(res, out)
cat("saved", out, length(res), "models\n")
