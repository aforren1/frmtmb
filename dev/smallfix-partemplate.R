# Review finding 8 says par_template() takes no prior. Its formula
# method does (par_template.default(prior =)); its fit method does not.
# Both, with a dirichlet() on class "b".
#   Rscript dev/smallfix-partemplate.R <lib;lib>
# Log: dev/smallfix-out/partemplate.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(x = rnorm(40))
d$y <- 1 + d$x + rnorm(40)
p <- set_prior("dirichlet(2)", class = "b")
say <- function(lab, expr) {
  r <- tryCatch({ force(expr); "accepted" }, error = function(e) {
    paste("REFUSED:", substr(gsub("\\s+", " ", conditionMessage(e)), 1, 150))
  })
  cat(lab, ":", r, "\n")
}
say("par_template(formula, prior =)",
    par_template(bf(y ~ x), data = d, family = gaussian(), prior = p))
fit <- frm(bf(y ~ x), family = gaussian(), data = d)
say("par_template(fit, prior =)", par_template(fit, prior = p))
