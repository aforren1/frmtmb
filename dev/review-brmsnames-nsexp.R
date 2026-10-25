## A namespaced special, base::exp(x), and other spellings the alias does
## not reach. Seed 4.  Rscript dev/review-brmsnames-nsexp.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
cat("arm", arm, dirname(system.file(package = "frmtmb")), "\n")
set.seed(4)
n <- 200
d <- data.frame(x = runif(n), z = rnorm(n))
d$y <- 1 + 2 * exp(d$x) + 0.5 * d$z + rnorm(n)
for (f in list(y ~ base::exp(x) + z, y ~ exp(exp(x / 3)) + z,
               y ~ z + exp(x)^2, y ~ log(exp(x)) + z)) {
  r <- tryCatch(rownames(fixef(q(frm(f, data = d)))),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  b <- tryCatch(colnames(q(brms::standata(f, data = d))$X),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-28s frm: %s\n%-28s brms: %s\n", deparse1(f),
              paste(r, collapse = " "), "", paste(b, collapse = " ")))
}
