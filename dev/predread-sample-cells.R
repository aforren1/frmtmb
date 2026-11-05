# frmtmb.sample's readers whose layouts the cs() change reaches, against brms
# 2.23.0 sampled on the same data (layout only; seed 4).
#   Rscript dev/predread-sample-cells.R <lib> > dev/predread-log/sample-cells.txt
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r16"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r16",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
quiet <- function(e) suppressWarnings(suppressMessages(e))
shape <- function(lab, e) {
  r <- tryCatch(quiet(e), error = function(err) {
    paste("ERROR:", substr(conditionMessage(err), 1, 150))
  })
  if (is.character(r) && length(r) == 1L) {
    cat(sprintf("%-46s %s\n", lab, r))
  } else {
    third <- if (length(dim(r)) == 3L) dimnames(r)[[3L]]
    cat(sprintf("%-46s dim %s third %s\n", lab, paste(dim(r), collapse = "x"),
                paste(third, collapse = " ")))
  }
}
set.seed(4)
d <- data.frame(x = stats::rnorm(150), z = stats::rnorm(150))
d$o <- cut(d$x - 0.5 * d$z + stats::rlogis(150), c(-Inf, -0.5, 0.5, Inf),
           labels = FALSE)
d$y <- d$x + stats::rnorm(150)
ds <- quiet(frm_sample(frmtmb::bf(o ~ x + cs(z)), data = d,
                       family = frmtmb::sratio(), chains = 1, iter = 300,
                       seed = 3, refresh = 0))
bs <- quiet(brm(brms::bf(o ~ x + cs(z)), data = d, family = brms::sratio(),
                chains = 1, iter = 300, seed = 3, refresh = 0))
nd <- data.frame(x = c(-1, 1), z = c(0.5, -1), y = 0)
for (a in list(
  list("fitted(scale = linear)", function(f) fitted(f, newdata = nd,
                                                    scale = "linear")),
  list("fitted(scale = linear, summary = FALSE)", function(f) {
    fitted(f, newdata = nd, scale = "linear", summary = FALSE)
  }),
  list("fitted(dpar = mu, scale = linear)", function(f) {
    fitted(f, newdata = nd, dpar = "mu", scale = "linear")
  }),
  list("posterior_linpred", function(f) posterior_linpred(f, newdata = nd)),
  list("posterior_linpred(incl_thres)", function(f) {
    posterior_linpred(f, newdata = nd, incl_thres = TRUE)
  }),
  list("loo_linpred", function(f) loo_linpred(f)),
  list("conditional_effects dpar = mu", function(f) {
    as.matrix(conditional_effects(f, "x", dpar = "mu",
                                  resolution = 3)[[1L]]["estimate__"])
  }))) {
  shape(paste("frmtmb", a[[1L]]), a[[2L]](ds))
  shape(paste("brms  ", a[[1L]]), a[[2L]](bs))
}
dm <- quiet(frm_sample(frmtmb::mvbf(frmtmb::bf(y ~ x), frmtmb::bf(o ~ x + cs(z),
                                                                  family = frmtmb::sratio()),
                                    rescor = FALSE),
                       data = d, chains = 1, iter = 200, seed = 3, refresh = 0))
bm <- quiet(brm(brms::mvbf(brms::bf(y ~ x), brms::bf(o ~ x + cs(z),
                                                     family = brms::sratio()),
                           rescor = FALSE),
                data = d, chains = 1, iter = 200, seed = 3, refresh = 0))
shape("frmtmb mv posterior_linpred", posterior_linpred(dm, newdata = nd))
shape("brms   mv posterior_linpred", posterior_linpred(bm, newdata = nd))
shape("frmtmb mv posterior_linpred resp = o",
      posterior_linpred(dm, newdata = nd, resp = "o"))
shape("brms   mv posterior_linpred resp = o",
      posterior_linpred(bm, newdata = nd, resp = "o"))
shape("frmtmb mv fitted linear", fitted(dm, newdata = nd, scale = "linear"))
shape("brms   mv fitted linear", fitted(bm, newdata = nd, scale = "linear"))
