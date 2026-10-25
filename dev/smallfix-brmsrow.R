# A brms prior object with a dirichlet() row on a class that is no
# simplex, passed to frmtmb (review finding 9): where it is refused and
# with what words, before and after.
#   Rscript dev/smallfix-brmsrow.R <lib;lib>
# Log: dev/smallfix-out/brmsrow-<arm>.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(1)
d <- data.frame(x = rnorm(40))
d$y <- 1 + d$x + rnorm(40)
bp <- brms::prior(dirichlet(2), class = b)
say <- function(lab, expr) {
  r <- tryCatch({ force(expr); "accepted" }, error = function(e) {
    paste("REFUSED:", gsub("\\s+", " ", conditionMessage(e)))
  })
  cat(lab, ":", r, "\n")
}
say("validate_prior(brms prior)", validate_prior(bp, bf(y ~ x), data = d))
say("frm(prior = brms prior)", frm(bf(y ~ x), family = gaussian(), data = d,
                                   prior = bp))
