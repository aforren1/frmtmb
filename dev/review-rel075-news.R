# Review: spot-check NEWS claims on rellib-r13, core attached alone.
.libPaths(c("C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
show <- function(label, expr) {
  cat("\n=====", label, "\n")
  r <- tryCatch(expr, error = function(e) {
    cat("ERROR [", paste(class(e), collapse = ","), "]:",
        conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(r)) print(r)
  invisible(r)
}
set.seed(1)
d <- data.frame(time = rexp(50), x = rnorm(50), c = 0)
show("bhaz without frmtmb.survival loaded",
     frm(time | bhaz(df = 6) ~ x, data = d, family = gaussian()))
cat("survival loaded after that call?", isNamespaceLoaded("frmtmb.survival"), "\n")
show("set_prior dirichlet on Intercept dpar sigma (accepted)",
     set_prior("dirichlet(2)", class = "Intercept", dpar = "sigma"))
show("validate_prior refuses it",
     validate_prior(set_prior("dirichlet(2)", class = "Intercept",
                              dpar = "sigma"), time ~ x, data = d))
show("validate_prior refuses dirichlet on class sigma",
     validate_prior(set_prior("dirichlet(2)", class = "sigma"), time ~ x,
                    data = d))
show("frmtmb_covstruct_info('car') loads frmtmb.spatial",
     {r <- frmtmb_covstruct_info("car"); c(r, loaded = isNamespaceLoaded("frmtmb.spatial"))})
show("frmtmb_covstruct_info('nope')", frmtmb_covstruct_info("nope"))
show("frmtmb_covstruct_info('us')", frmtmb_covstruct_info("us"))
fit <- frm(time ~ x, data = d)
cv <- frm_curve(fit, newdata = data.frame(x = seq(-1, 1, length.out = 5)))
show("print(curve, digits = 2)", {print(cv, digits = 2); "ok"})
show("cox without survival",
     frm(time | cens(c) ~ x, data = d, family = "cox"))
