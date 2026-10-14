.libPaths(c("C:/Users/adf44/source/r/wt-famsimplex-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
dd <- data.frame(g = factor(rep(c("a", "b"), each = 5)))
dd$Y <- matrix(1/3, 10, 3, dimnames = list(NULL, c("a", "b", "c")))
str(unclass(par_template(Y ~ 0 + g, dd, family = dirichlet())))
str(unclass(par_template(Y ~ 0 + g, dd, family = logistic_normal())))
str(unclass(par_template(Y ~ 0 + g, dd, family = brmsfamily("dirichlet2"))))
s <- frm_simulate(Y ~ 0 + g, dd, family = dirichlet(),
                  newparams = list(b_mub_ga = 0.5, b_mub_gb = -0.4,
                                   b_muc_ga = -0.3, b_muc_gb = 0.6,
                                   phi = 8), nsim = 1, seed = 1)
str(s)
