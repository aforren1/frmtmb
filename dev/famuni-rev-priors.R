# Reviewer probe: brms 2.23.0 default_prior() against frmtmb's sampling
# route for the three families, and the scalar-xi note.
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb", format(packageVersion("frmtmb")), "from", find.package("frmtmb"),
    "; sample from", find.package("frmtmb.sample"), "; brms",
    format(packageVersion("brms")), "\n")
set.seed(20261017)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$yf <- exp(1 + 0.3 * d$x) / gamma(1 - 1 / 4) * rexp(n)^(-1 / 4)
d$yg <- 1 + 0.5 * d$x + 0.8 * (rexp(n)^(-0.15) - 1) / 0.15
d$yd <- pmax(ceiling((log1p(-runif(n)) / log(plogis(1 + 0.3 * d$x)))^(1 / 1.2) - 1), 0)
cases <- list(
  list("frechet y ~ x", "yf ~ x", "frechet", NULL, list()),
  list("frechet nu ~ z", "yf ~ x", "frechet", "nu ~ z", list()),
  list("GEV sigma ~ z, xi ~ z", "yg ~ x", "gen_extreme_value", c("sigma ~ z", "xi ~ z"), list()),
  list("GEV scalar xi", "yg ~ x", "gen_extreme_value", "sigma ~ z", list()),
  list("GEV xi ~ z, identity", "yg ~ x", "gen_extreme_value", "xi ~ z", list(link_xi = "identity")),
  list("GEV xi ~ 1 (intercept only)", "yg ~ x", "gen_extreme_value", "xi ~ 1", list()),
  list("dw y ~ x", "yd ~ x", "discrete_weibull", NULL, list()),
  list("dw shape ~ z", "yd ~ x", "discrete_weibull", "shape ~ z", list()))
for (cs in cases) {
  cat("\n=====", cs[[1]], "\n")
  bfo <- do.call(brms::bf, c(list(as.formula(cs[[2]])),
                             lapply(cs[[4]], as.formula)))
  bfam <- do.call(brms::brmsfamily, c(list(cs[[3]]), cs[[5]]))
  bp <- suppressWarnings(brms::default_prior(bfo, data = d, family = bfam))
  bp <- as.data.frame(bp)
  bp <- bp[nzchar(bp$prior), c("prior", "class", "coef", "dpar", "lb")]
  cat("brms:\n"); print(bp, row.names = FALSE)
  ffo <- do.call(frmtmb::bf, c(list(as.formula(cs[[2]])),
                               lapply(cs[[4]], as.formula)))
  ffam <- do.call(get(cs[[3]], asNamespace("frmtmb")), cs[[5]])
  fp <- tryCatch(frmtmb::default_prior(ffo, data = d, family = ffam,
                                       route = "sample"),
                 error = function(e) conditionMessage(e))
  cat("frmtmb route = 'sample':\n")
  if (is.character(fp)) cat(fp, "\n") else {
    fp <- as.data.frame(fp)
    keep <- intersect(c("prior", "class", "coef", "dpar", "source"), names(fp))
    print(fp[nzchar(fp$prior), keep], row.names = FALSE)
  }
  fit <- frm(ffo, family = ffam, data = d)
  notes <- frmtmb.sample:::default_prior_notes(fit)
  cat("notes:\n"); cat(paste0("  - ", notes), sep = "\n")
}
