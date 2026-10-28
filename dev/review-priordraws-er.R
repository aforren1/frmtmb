# Reviewer: Evid.Ratio of frmtmb's hypothesis() against brms's own
# evidence_ratio() and prior_draws() sequence on the same draws: brms
# reads prior_draws(x, parsH) for EVERY hypothesis, directional or not,
# one sample() per variable, after set.seed(seed).
#   Rscript dev/review-priordraws-er.R
.libPaths(c("C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-stancache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
options(mc.cores = 1)
set.seed(1)
d <- data.frame(x = rnorm(80), z = rnorm(80), g = factor(rep(1:8, 10)))
d$y <- 1 + 0.5 * d$x + 0.2 * d$z + rnorm(8, 0, 0.7)[d$g] + rnorm(80)
yes <- suppressWarnings(suppressMessages(frm_sample(
  bf(y ~ x + z + (1 | g)), data = d, family = gaussian(),
  prior = set_prior("normal(0, 1)", class = "b"), sample_prior = "yes",
  chains = 2, iter = 1000, seed = 21, refresh = 0)))
dr <- as.data.frame(yes)
find_vars <- utils::getFromNamespace("find_vars", "brms")
evr <- utils::getFromNamespace("evidence_ratio", "brms")
ref <- function(hs, class, seed) {
  set.seed(seed)
  sapply(hs, function(h) {
    h <- gsub("[ \t\r\n]", "", h)
    sign <- regmatches(h, regexpr("=|<|>", h))
    lr <- strsplit(h, "=|<|>")[[1]]
    e <- paste0("(", lr[1], ")", if (lr[2] != "0") paste0("-(", lr[2], ")"))
    v <- find_vars(e)
    pars <- paste0(class, v)
    post <- eval(parse(text = e), setNames(as.list(dr[pars]), v))
    pd <- prior_draws(yes, variable = pars)
    pri <- if (ncol(pd) == length(v)) {
      eval(parse(text = e), setNames(as.list(pd), v))
    }
    evr(post, cut = 0, wsign = switch(sign, "=" = "equal", "<" = "less",
                                      ">" = "greater"),
        prior_samples = pri)
  })
}
hs <- c("x = 0", "x > 0", "x < 0.3", "x - z = 0", "x + z = 0.5", "z = x")
fr <- hypothesis(yes, hs, seed = 3)$hypothesis$Evid.Ratio
br <- ref(hs, "b_", 3)
for (k in seq_along(hs)) cat(sprintf("  %-12s frmtmb %.12g  brms %.12g  identical %s\n",
                                     hs[k], fr[k], br[k], identical(unname(fr[k]), unname(br[k]))))
hs2 <- c("g__Intercept = 0.5", "g__Intercept > 1")
fr2 <- tryCatch(hypothesis(yes, hs2, class = "sd", seed = 4)$hypothesis$Evid.Ratio,
                error = function(e) conditionMessage(e))
br2 <- ref(hs2, "sd_", 4)
cat("class = sd: frmtmb", fr2, " brms", br2, "\n")
hs3 <- c("sigma = 1", "sigma^2 = 1")
fr3 <- hypothesis(yes, hs3, class = "", seed = 5)$hypothesis$Evid.Ratio
br3 <- ref(hs3, "", 5)
cat("class = '': frmtmb", fr3, " brms", br3, "\n")
