# What brms 2.23.0 does with cs() on cumulative(), on
# hurdle_cumulative() and on a cumulative() component of an ordinal
# mixture: the warnings it gives, the Stan likelihood it writes and the
# thresholds it reads, run rather than recalled (dev/round-20261005.md,
# small change a).
#
#   Rscript dev/rel068-cs-brms.R > dev/rel068-log/cs-brms.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
d$y <- sample(1:4, n, TRUE)
d$yh <- sample(0:4, n, TRUE)
run <- function(tag, f, fam) {
  cat("\n==", tag, "\n")
  w <- character()
  code <- tryCatch(withCallingHandlers(
    as.character(stancode(f, data = d, family = fam)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage")),
    error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("warnings:", length(w), "\n")
  for (x in w) cat("  W:", x, "\n")
  if (startsWith(code, "ERROR")) {
    cat(code, "\n")
    return(invisible())
  }
  l <- strsplit(code, "\n")[[1]]
  keep <- grep("lpmf|lupmf|ordered|Intercept|bcs|mucs", l, value = TRUE)
  cat(paste0("  ", trimws(keep)), sep = "\n")
}
run("cumulative logit, y ~ cs(x)", bf(y ~ cs(x)), cumulative())
run("cumulative probit, y ~ z + cs(x), disc ~ 0 + z",
    bf(y ~ z + cs(x), disc ~ 0 + z), cumulative("probit"))
run("cumulative equidistant, y ~ cs(x)", bf(y ~ cs(x)),
    cumulative(threshold = "equidistant"))
run("hurdle_cumulative probit, yh ~ cs(x)", bf(yh ~ cs(x)),
    hurdle_cumulative("probit"))
run("mixture(cumulative, cumulative), y ~ cs(x)", bf(y ~ cs(x)),
    mixture(cumulative(), cumulative()))
run("mixture(cumulative, sratio), y ~ cs(x)", bf(y ~ cs(x)),
    mixture(cumulative(), sratio()))
run("sratio, y ~ cs(x) (control: no warning expected)", bf(y ~ cs(x)),
    sratio())
run("cumulative, y ~ x (control: no cs)", bf(y ~ x), cumulative())
cat("\nbrms's cumulative_logit_lpmf body:\n")
code <- as.character(stancode(bf(y ~ x), data = d, family = cumulative()))
l <- strsplit(code, "\n")[[1]]
s <- grep("real cumulative_logit_lpmf", l)
if (length(s)) cat(l[s:(s + 14)], sep = "\n")
