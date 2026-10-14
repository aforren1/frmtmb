# Lane dblbar: the post-fit paths a `||` term reaches (update(), refit(),
# frmtmb.sample's draws and priors, as_tmbstan(), posterior_samples()
# names), on (1 + f || g) and its explicit gr(cor = FALSE) spelling.
# Data seed 20261013 (dev/dblbar-data.R); sampler seed 1.
# Usage: Rscript dev/dblbar-smoke.R base|lane
# Output: dev/dblbar-log/smoke-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-dblbar/dev/stan-cache")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
dd <- dblbar_data()
sink(file.path(wt, sprintf("dev/dblbar-log/smoke-%s.txt", arm)),
     split = TRUE)
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
try_show <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:",
                                                conditionMessage(e)))
  cat("--", label, ":", if (is.character(r)) r else "", "\n")
  if (!is.character(r)) print(r)
  invisible(r)
}
fit <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
try_show("variables", variables(fit))
u1 <- update(fit, formula. = ~ . + x)
try_show("update adds x: sd names",
         grep("^sd_|^cor_", variables(u1), value = TRUE))
u2 <- update(fit, newdata = dd[1:300, ])
try_show("update newdata: blocks", vapply(u2$frame$re_blocks, function(b) {
  paste0(b$covstruct, "[", paste(b$cnms, collapse = ","), "]")
}, ""))
r1 <- tryCatch(refit(fit, data = dd[1:300, ]), error = function(e) e)
try_show("refit identical to update", if (inherits(r1, "error")) {
  conditionMessage(r1)
} else identical(r1$estimates, u2$estimates))
try_show("default_prior route sample: sd/cor rows", {
  p <- as.data.frame(default_prior(bf(y ~ f + (1 + f || g)), data = dd,
                                   family = gaussian(), route = "sample"))
  p[p$class %in% c("sd", "cor"), c("prior", "class", "coef", "group")]
})
try_show("prior on coef fb is accepted", {
  f2 <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd,
            prior = set_prior("normal(0, 0.1)", class = "sd", coef = "fb",
                              group = "g"))
  c(sd_fb_flat = unname(VarCorr(fit)$g$sd["fb", 1]),
    sd_fb_prior = unname(VarCorr(f2)$g$sd["fb", 1]))
})
set.seed(1)
dr <- tryCatch(frm_sample(fit, chains = 2, iter = 600, seed = 1,
                          refresh = 0),
               error = function(e) e)
if (inherits(dr, "error")) {
  cat("frm_sample ERROR:", conditionMessage(dr), "\n")
} else {
  try_show("draws variables (sd, cor)",
           grep("^sd_|^cor_", variables(dr), value = TRUE))
  try_show("draws VarCorr g elements", names(VarCorr(dr)$g))
  try_show("draws summary random rows", rownames(summary(dr)$random$g))
  try_show("posterior_samples sd names",
           grep("^sd_|^cor_", names(suppressWarnings(
             posterior_samples(dr))), value = TRUE))
  try_show("prior_summary sd/cor rows", {
    ps <- as.data.frame(prior_summary(dr))
    ps[ps$class %in% c("sd", "cor"), c("prior", "class", "coef", "group")]
  })
}
ts <- tryCatch(as_tmbstan(fit, chains = 1, iter = 300, seed = 1,
                          refresh = 0),
               error = function(e) e)
try_show("as_tmbstan", if (inherits(ts, "error")) conditionMessage(ts) else
  class(ts))
sink()
