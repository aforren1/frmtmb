# brms 2.23.0 on the item-3 construction: predict(newdata = ) with the
# response NA in rows 5 to 8 of group 1, under ar(cov = FALSE). Same
# data as dev/formrobust-repro3.R (seed 31). Stan seed 3, 1 chain, 2000
# iterations.
.libPaths(c("C:/Users/adf44/source/r/wt-formrobust-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
set.seed(31)
G <- 30; Tn <- 8
d <- expand.grid(t = 1:Tn, g = factor(1:G))
d$x <- rnorm(nrow(d))
e <- as.vector(apply(matrix(rnorm(G * Tn), Tn, G), 2, function(z) {
  as.vector(stats::filter(z, 0.6, "recursive"))
}))
d$y <- 1 + 0.5 * d$x + e
nd <- d[d$g %in% c("1", "2"), ]
nd$y[nd$g == "1" & nd$t >= 5] <- NA
fb <- brms::brm(y ~ x + ar(t, g, p = 1), data = d, chains = 1, iter = 2000,
                seed = 3, refresh = 0)
set.seed(1)
pb <- predict(fb, newdata = nd)
fe <- fitted(fb, newdata = nd)
i5 <- which(nd$g == "1" & nd$t == 5)
i6 <- i5 + 1L
suppressMessages(library(frmtmb))
ff <- frm(bf(y ~ x + ar(t, g, p = 1)), data = d)
set.seed(1)
pf <- predict(ff, newdata = nd, ndraws = 4000)
ffit <- fitted(ff, newdata = nd)
cat("brms   predict Est.Error rows 5..8 of g1:",
    format(pb[i5:(i5 + 3), "Est.Error"], digits = 5), "\n")
cat("frmtmb predict Est.Error rows 5..8 of g1:",
    format(pf[i5:(i5 + 3), "Est.Error"], digits = 5), "\n")
cat("brms   fitted Estimate rows 5..8 of g1:",
    format(fe[i5:(i5 + 3), "Estimate"], digits = 5), "\n")
cat("frmtmb fitted Estimate rows 5..8 of g1:",
    format(ffit[i5:(i5 + 3), "Estimate"], digits = 5), "\n")
cat("brms   fitted Est.Error rows 5..8 of g1:",
    format(fe[i5:(i5 + 3), "Est.Error"], digits = 5), "\n")
cat("frmtmb fitted Est.Error rows 5..8 of g1:",
    format(ffit[i5:(i5 + 3), "Est.Error"], digits = 5), "\n")
cat("brms ar:", format(brms::fixef(fb)[, 1], digits = 5),
    format(summary(fb)$cor_pars[, 1], digits = 5), "\n")
