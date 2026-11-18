# Reviewer of lane ordstarts, build round. Probe 7 (question 2): K = 3 with
# the weights' reference moved off the last component (formulas on theta2
# and theta3): does mixture_ord_swap_starts() still build its starts?
#   Rscript ordstarts-build-review-p7.R > ordstarts-build-review-p7.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
set.seed(3)
n <- 900
d <- data.frame(x = rnorm(n), z = rnorm(n))
cls <- sample(1:3, n, TRUE)
lat <- c(1.5, -0.8, 0)[cls] * d$x + c(1, -1, 0)[cls] + rlogis(n)
d$y <- 1L + (lat > -1) + (lat > 0) + (lat > 1)
fam <- mixture(cumulative(), sratio(), acat())
q <- function(e) suppressWarnings(suppressMessages(e))
for (form in list(bf(y ~ x), bf(y ~ x, theta1 ~ 1, theta2 ~ 1),
                  bf(y ~ x, theta2 ~ 1, theta3 ~ 1))) {
  f <- q(frm(form, family = fam, data = d,
             control = frmtmb_control(starts = 0)))
  nm <- ns$outer_par_map(f)$names
  sw <- tryCatch(ns$mixture_ord_swap_starts(f, "y"), error = function(e) e)
  dd <- q(frm(form, family = fam, data = d))
  ms <- dd$opt$mixture_starts
  cat(sprintf("%-40s weights %s | swap starts built: %s | default record swaps %s\n",
              paste(names(form$pforms), collapse = ","),
              paste(grep("^theta", nm, value = TRUE), collapse = ","),
              if (inherits(sw, "error")) paste("ERROR", conditionMessage(sw))
              else length(sw),
              if (is.null(ms)) "NULL" else ms[["swaps"]]))
}
cat("DONE p7\n")
