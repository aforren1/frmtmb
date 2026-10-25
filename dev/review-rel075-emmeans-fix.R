# Review: the proposed fix for recover_data.frmtmb_fit(), tried in this
# process only (the namespace copy is patched in memory; no file changes):
# drop the head check, so emm_basis.frmtmb_fit()'s head check raises the
# classed refusal outside emmeans's try(), as on 0.74.0.
.libPaths(c("C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
f <- get("recover_data.frmtmb_fit", ns)
b <- body(f)
stopifnot(identical(b[[2]], quote(require_fit_packages(object))))
body(f) <- b[-2]
unlockBinding("recover_data.frmtmb_fit", ns); assign("recover_data.frmtmb_fit", f, envir = ns)
registerS3method("recover_data", "frmtmb_fit", f, envir = asNamespace("emmeans"))
set.seed(1)
d <- data.frame(g = factor(rep(1:20, each = 5)), x = stats::rnorm(100),
                f = factor(rep(c("a", "b"), 50)))
d$y <- 1 + d$x + (d$f == "b") + stats::rnorm(20)[d$g] + stats::rnorm(100, 0, 0.5)
fit <- frm(y ~ x + f + (1 | g), data = d)
print(emmeans::emmeans(fit, ~ f))
fam <- fit$spec$responses[[1L]]$family
attr(fam, "frmtmb_package") <- "frmtmbabsenttestpkga"
fit$spec$responses[[1L]]$family <- fam
fit$frame$spec <- fit$spec
fit$cache <- new.env()
for (nm in c("ref_grid", "emmeans")) {
  r <- tryCatch({ if (nm == "ref_grid") emmeans::ref_grid(fit) else
    emmeans::emmeans(fit, ~ f); "ANSWERED" },
    error = function(e) paste0("ERROR [", class(e)[1], "]: ", conditionMessage(e)))
  cat(nm, ":", r, "\n")
}
