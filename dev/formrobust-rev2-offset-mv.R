# Re-check: a multivariate fit of one family where only one response
# has an offset, which goes through the new grid-route switch; brms at
# fixed parameters. Data seed 21, as formrobust-rev2-offset.R.
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev2-offset.R")
eval(parse(text = src[1:grep("^mx <- mean", src)]))
cat("\n== E. mv poisson, offset on yc only ==\n")
fE <- frm(bf(yc ~ x + f + offset(log(time))) + bf(yc2 ~ x + f), data = d,
          family = poisson())
bE <- tr("brms fixed E", brms_fixed_fit(
  brms::bf(yc ~ x + f + offset(log(time))) + brms::bf(yc2 ~ x + f) +
    brms::set_rescor(FALSE), brms::brmsfamily("poisson"), d, fE,
  ndraws = 4))
be <- fixef(fE)[, "Estimate"]
specs <- list(list("both", list()),
              list("both at time c(1,2)", list(at = list(time = c(1, 2)))),
              list("both | time at c(1,2)", list(at = list(time = c(1, 2)))),
              list("both epred", list(epred = TRUE)),
              list("resp yc2", list(resp = "yc2")),
              list("both by h? none", list()))
for (sp in specs[1:5]) {
  frm_spec <- if (grepl("[|] time", sp[[1]])) ~ f | time else ~ f
  ef <- tr(paste("frmtmb E", sp[[1]]),
           summary(do.call(emmeans, c(list(fE, frm_spec), sp[[2]]))))
  eb <- if (!is.null(bE)) tr(paste("brms E", sp[[1]]),
                             summary(do.call(emmeans,
                                             c(list(bE, frm_spec), sp[[2]]))))
  cmp(paste("mvE", sp[[1]]), ef, eb)
  if (!is.null(ef)) print(as.data.frame(ef)[, 1:3])
  if (!is.null(eb)) print(as.data.frame(eb)[, 1:3])
}
lin1 <- be[["yc_Intercept"]] + be[["yc_x"]] * mx + c(0, be[["yc_fb"]])
lin2 <- be[["yc2_Intercept"]] + be[["yc2_x"]] * mx + c(0, be[["yc2_fb"]])
cat("hand yc: lin + lmt", format(lin1 + lmt, digits = 10), "; yc2: lin",
    format(lin2, digits = 10), "\n")
cat("hand at time c(1,2): yc", format(lin1 + mean(log(c(1, 2))), digits = 10),
    "\n")
