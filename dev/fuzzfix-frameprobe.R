.libPaths(c("C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
d <- readRDS("dev/fuzzfix-out/data-gam.rds")
fr <- frm(bf(y | weights(w) ~ 1 + x + mo(xo) * z + ar1(tt + 0 | g), shape ~ 1 + z),
          family = Gamma(link = "log"), data = d, dry_run = "frame")
str(names(fr))
for (k in names(fr$linpreds)) {
  lp <- fr$linpreds[[k]]
  cat(k, "resp", lp$resp, "dpar", lp$dpar, "Z", if (!is.null(lp$Z)) dim(lp$Z), class(lp$Z), "const", deparse(lp$constant), "nl", !is.null(lp$nl_body), "\n")
}
str(fr$aterm_values, max.level = 2)
str(fr$re_blocks[[1]], max.level = 1)
cat("n_obs", fr$n_obs, "\n")
sp <- frm(bf(y | weights(w) ~ 1 + x + ar1(tt + 0 | g), shape ~ 1 + z),
          family = Gamma(link = "log"), data = d, dry_run = "spec")
str(sp$responses[[1]]$family[c("family", "dpars", "primary_dpars", "type")])
str(names(sp$responses[[1]]))
