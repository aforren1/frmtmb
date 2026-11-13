# Reviewer of lane priortab: sigma links softplus/squareplus, `+` identity.
.libPaths(c("C:/Users/adf44/source/r/wt-priortab-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({loadNamespace("brms"); loadNamespace("frmtmb")})
F <- asNamespace("frmtmb"); B <- asNamespace("brms")
cat("identical + methods:", identical(F[["+.frmtmb_prior_rows"]],
                                      F[["+.frmtmb_priorlist"]]), "\n")
set.seed(1)
d <- data.frame(z = rnorm(50)); d$y <- d$z + rnorm(50)
for (lk in c("softplus", "squareplus", "identity", "log")) {
  b <- tryCatch(as.data.frame(B$default_prior(y ~ z, data = d,
         family = B$brmsfamily("gaussian", link_sigma = lk))),
         error = function(e) conditionMessage(e))
  f <- tryCatch(as.data.frame(F$default_prior(y ~ z, data = d,
         family = F$brmsfamily("gaussian", link_sigma = lk))),
         error = function(e) conditionMessage(e))
  cat(sprintf("link_sigma = %-10s brms sigma lb: %-4s frmtmb sigma lb: %s\n", lk,
              if (is.data.frame(b)) b$lb[b$class == "sigma"] else b,
              if (is.data.frame(f)) f$lb[f$class == "sigma"] else substr(f, 1, 80)))
}
for (lk in c("softplus", "squareplus", "identity", "log")) {
  b <- tryCatch(as.data.frame(B$default_prior(y ~ z, data = d,
         family = B$brmsfamily("student", link_nu = if (lk == "log") "logm1" else lk))),
         error = function(e) conditionMessage(e))
  f <- tryCatch(as.data.frame(F$default_prior(y ~ z, data = d,
         family = F$brmsfamily("student", link_nu = if (lk == "log") "logm1" else lk))),
         error = function(e) conditionMessage(e))
  cat(sprintf("link_nu = %-10s brms nu lb: %-4s frmtmb nu lb: %s\n", lk,
              if (is.data.frame(b)) b$lb[b$class == "nu"] else substr(b, 1, 60),
              if (is.data.frame(f)) f$lb[f$class == "nu"] else substr(f, 1, 80)))
}
