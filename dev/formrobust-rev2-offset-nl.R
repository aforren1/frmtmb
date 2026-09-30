# Re-check: nonlinear model with an offset in nlpar a; brms emmeans at
# fixed parameters for the body, nlpar = "a" and epred. Data seed 21.
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev2-offset.R")
eval(parse(text = src[1:grep("^mx <- mean", src)]))
d$y2n <- rnorm(n, 1 + 0.5 * d$x + d$z, 0.5)
fN <- frm(bf(y2n ~ a + b * x, a ~ 1 + offset(z), b ~ 1, nl = TRUE), data = d)
bN <- tr("brms fixed N", brms_fixed_fit(
  brms::bf(y2n ~ a + b * x, a ~ 1 + offset(z), b ~ 1, nl = TRUE),
  brms::brmsfamily("gaussian"), d, fN, ndraws = 4))
for (sp in list(list("body", list()), list("nlpar a", list(nlpar = "a")),
                list("nlpar b", list(nlpar = "b")),
                list("epred", list(epred = TRUE)),
                list("nlpar a at z = 0", list(nlpar = "a", at = list(z = 0))))) {
  ef <- tr(paste("frmtmb", sp[[1]]), summary(do.call(emmeans,
                                                    c(list(fN, ~ 1), sp[[2]]))))
  eb <- if (!is.null(bN)) tr(paste("brms", sp[[1]]),
                             summary(do.call(emmeans, c(list(bN, ~ 1), sp[[2]]))))
  cmp(paste("nl", sp[[1]]), ef, eb)
}
