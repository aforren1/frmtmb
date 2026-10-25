# Review: what a user sees from emmeans on a fit with a part this
# session cannot use. recover_data.frmtmb_fit() calls
# require_fit_packages() OUTSIDE its tryCatch(), and its own comment says
# emmeans calls it inside try() and replaces whatever it raises.
# The fit is built as test-extension-home-mixed.R builds it: a family
# said to be owned by a package no session has.
.libPaths(c("C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(g = factor(rep(1:20, each = 5)), x = stats::rnorm(100),
                f = factor(rep(c("a", "b"), 50)))
d$y <- 1 + d$x + (d$f == "b") + stats::rnorm(20)[d$g] +
  stats::rnorm(100, 0, 0.5)
fit <- frm(y ~ x + f + (1 | g), data = d)
cat("intact fit, emmeans(fit, ~ f):\n")
print(emmeans::emmeans(fit, ~ f))
fam <- fit$spec$responses[[1L]]$family
attr(fam, "frmtmb_package") <- "frmtmbabsenttestpkga"
fit$spec$responses[[1L]]$family <- fam
fit$frame$spec <- fit$spec
fit$cache <- new.env()
show <- function(label, expr) {
  r <- tryCatch({ force(expr); "ANSWERED" }, error = function(e)
    paste0("ERROR [", paste(class(e), collapse = ","), "]: ",
           conditionMessage(e)))
  cat("\n==", label, "\n", r, "\n")
}
show("summary(fit)", summary(fit))
show("emmeans::ref_grid(fit)", emmeans::ref_grid(fit))
show("emmeans::emmeans(fit, ~ f)", emmeans::emmeans(fit, ~ f))
show("recover_data() called directly",
     emmeans::recover_data(fit))
show("marginaleffects::avg_slopes(fit)",
     marginaleffects::avg_slopes(fit))
show("marginaleffects::predictions(fit)",
     marginaleffects::predictions(fit))
