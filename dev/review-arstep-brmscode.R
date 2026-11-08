# Review of lane arstep: does brms 2.23.0 bound the cov = FALSE ar
# coefficient? stancode() only, no compile.
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-lib.R")
rv_libs("lane")
suppressPackageStartupMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), dirname(find.package("brms")), "\n")
set.seed(1)
n <- 40
d <- data.frame(x = sort(runif(n, 0, 10)), t = seq_len(n))
d$y <- 2 * (d$x > 5) - 1 + rnorm(n, 0, 0.05)
for (cv in c(FALSE, TRUE)) {
  sc <- brms::stancode(y ~ gp(x) + ar(time = t, cov = cv), data = d)
  l <- grep("ar;|ar\\[|Kar\\]", strsplit(sc, "\n")[[1]], value = TRUE)
  cat("BRMS cov =", cv, "\n"); cat(paste(" ", l), sep = "\n")
  pr <- brms::default_prior(y ~ gp(x) + ar(time = t, cov = cv), data = d)
  print(pr[pr$class == "ar", c("prior", "class", "lb", "ub")])
}
cat("DONE\n")
