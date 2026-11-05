# What brms 2.23.0 returns from every reader of an ordinal predictor with a
# cs() term: layout only, on short sampled fits (seed 13).
#   Rscript dev/predread-cs-layout-brms.R > dev/predread-log/cs-layout-brms.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
set.seed(13)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- cut(d$x + 0.5 * d$z + rlogis(n), c(-Inf, -1, 0, 1.5, Inf),
           labels = FALSE)
nd <- data.frame(x = c(-1, 0, 1), z = c(0.2, 0.2, 0.2))
tr <- function(lab, e) {
  cat("\n---", lab, "\n")
  r <- tryCatch(e, error = function(err) {
    cat("ERROR:", conditionMessage(err), "\n")
    NULL
  })
  if (is.array(r) || is.matrix(r)) {
    cat("dim:", dim(r), "\n")
    str(dimnames(r))
    if (length(dim(r)) == 3L) print(round(r[1, , , drop = TRUE], 4)) else
      print(round(r[1, ], 4))
  } else if (!is.null(r)) {
    print(utils::head(as.data.frame(r), 8))
  }
  invisible(r)
}
quiet <- function(e) suppressWarnings(suppressMessages(e))
f1 <- quiet(brm(bf(y ~ x + cs(z)), data = d, family = sratio(), chains = 1,
                iter = 200, refresh = 0, seed = 1))
tr("sratio cs: posterior_linpred", posterior_linpred(f1, newdata = nd,
                                                     re_formula = NA))
tr("sratio cs: posterior_linpred(dpar = 'mu')",
   posterior_linpred(f1, newdata = nd, re_formula = NA, dpar = "mu"))
tr("sratio cs: posterior_linpred(transform = TRUE)",
   posterior_linpred(f1, newdata = nd, re_formula = NA, transform = TRUE))
tr("sratio cs: posterior_linpred(incl_thres = TRUE)",
   posterior_linpred(f1, newdata = nd, re_formula = NA, incl_thres = TRUE))
tr("sratio cs: posterior_epred(dpar = 'mu')",
   posterior_epred(f1, newdata = nd, re_formula = NA, dpar = "mu"))
tr("sratio cs: fitted(dpar = 'mu', scale = 'linear')",
   fitted(f1, newdata = nd, re_formula = NA, dpar = "mu", scale = "linear"))
tr("sratio cs: fitted(dpar = 'mu')",
   fitted(f1, newdata = nd, re_formula = NA, dpar = "mu"))
tr("sratio cs: emmeans ~ x", as.data.frame(emmeans::emmeans(
  f1, ~ x, at = list(x = c(-1, 1)))))
tr("sratio cs: emmeans ~ x | rep.meas", as.data.frame(emmeans::emmeans(
  f1, ~ x | rep.meas, at = list(x = c(-1, 1)))))
tr("sratio cs: conditional_effects default",
   conditional_effects(f1, "x")[[1]][, c("x", "estimate__", "se__",
                                         "lower__", "upper__")])
tr("sratio cs: conditional_effects(dpar = 'mu')",
   conditional_effects(f1, "x", dpar = "mu")[[1]][, c("x", "estimate__",
                                                      "se__")])
tr("sratio cs: conditional_effects(method = 'posterior_linpred')",
   conditional_effects(f1, "x", method = "posterior_linpred")[[1]][
     , c("x", "estimate__", "se__")])
tr("sratio cs: marginaleffects predictions type = 'link'",
   marginaleffects::predictions(f1, newdata = nd, type = "link"))
tr("sratio cs: marginaleffects predictions default",
   marginaleffects::predictions(f1, newdata = nd))

# a mixture and a hurdle
f2 <- quiet(brm(bf(y ~ 1, mu1 ~ x + cs(z), mu2 ~ x), data = d,
                family = mixture(cumulative(), cumulative()), chains = 1,
                iter = 200, refresh = 0, seed = 1, init = 0))
tr("mixture cs: posterior_linpred",
   posterior_linpred(f2, newdata = nd, re_formula = NA))
tr("mixture cs: posterior_linpred(dpar = 'mu1')",
   posterior_linpred(f2, newdata = nd, re_formula = NA, dpar = "mu1"))
tr("mixture cs: posterior_linpred(dpar = 'mu2')",
   posterior_linpred(f2, newdata = nd, re_formula = NA, dpar = "mu2"))
tr("mixture cs: fitted(scale = 'linear')",
   fitted(f2, newdata = nd, re_formula = NA, scale = "linear"))

f3 <- quiet(brm(bf(y | thres(gr = g) ~ x + cs(z)),
                data = transform(d, g = rep(c("a", "b"), 100)),
                family = cumulative(), chains = 1, iter = 200, refresh = 0,
                seed = 1, init = 0))
tr("thres(gr) cs: posterior_linpred",
   posterior_linpred(f3, newdata = transform(nd, g = c("a", "b", "a")),
                     re_formula = NA))
d4 <- d
d4$y <- ifelse(d4$y == 1L, 0L, d4$y - 1L)
f4 <- quiet(brm(bf(y ~ x + cs(z)), data = d4, family = hurdle_cumulative(),
                chains = 1, iter = 200, refresh = 0, seed = 1, init = 0))
tr("hurdle_cumulative cs: posterior_linpred",
   posterior_linpred(f4, newdata = nd, re_formula = NA))
