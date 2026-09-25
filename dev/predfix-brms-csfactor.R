# brms 2.23.0 on cs(factor), the wt-predfix reviewer's r2-cs2.R data,
# seed 405: does brms build dummy columns, treat the factor as its
# integer codes, or refuse? Parsing only (make_standata), no sampling.
#   Rscript dev/predfix-brms-csfactor.R > dev/predfix-log/brms-csfactor.txt
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(405)
n <- 500
x <- rnorm(n)
fc <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(a = -1, b = 0, c = 1.5)[as.character(fc)]
p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n)
d <- data.frame(x, fc, yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
r <- tryCatch({
  sd <- make_standata(bf(yo ~ cs(fc)), family = sratio(), data = d)
  cat("standata names:", paste(names(sd), collapse = ", "), "\n")
  cat("Kcs =", sd$Kcs, "; Xcs columns:", paste(colnames(sd$Xcs), collapse = ", "), "\n")
  print(head(sd$Xcs))
  "built"
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("make_standata:", r, "\n")
cat(tryCatch(paste(capture.output(get_prior(bf(yo ~ cs(fc)), family = sratio(),
                                            data = d)), collapse = "\n"),
             error = function(e) paste("get_prior ERROR:", conditionMessage(e))), "\n")
