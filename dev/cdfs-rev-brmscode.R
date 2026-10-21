# Reviewer: brms 2.23.0's Stan code for the censored and truncated rows
# of the discrete families, and its R-side truncation checks.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
d <- data.frame(y = rpois(20, 3), x = rnorm(20), cc = rep(c(-1, 0, 1, 2), 5),
                lb = 0, ub = 30, tr = 10)
d$y2 <- d$y + 2
show <- function(f, fam) {
  sc <- stancode(f, data = d, family = fam)
  ln <- strsplit(sc, "\n")[[1]]
  keep <- grep("cens|lccdf|lcdf|log_diff_exp|T\\[|lb|ub|Ycens|rcens",
               ln, value = TRUE)
  keep <- keep[!grepl("^\\s*//", keep)]
  cat("----", deparse(f), "|", fam$family, "\n")
  cat(unique(trimws(keep)), sep = "\n")
}
fams <- list(poisson(), negbinomial(), geometric(), binomial(),
             beta_binomial(), brmsfamily("com_poisson"), zero_inflated_poisson(),
             hurdle_poisson(), zero_inflated_negbinomial(),
             hurdle_negbinomial(), zero_inflated_binomial(),
             zero_inflated_beta_binomial())
for (fm in fams) {
  f1 <- if (fm$family %in% c("binomial", "beta_binomial",
                             "zero_inflated_binomial",
                             "zero_inflated_beta_binomial")) {
    y | trials(tr) + cens(cc, y2) ~ x
  } else y | cens(cc, y2) ~ x
  show(f1, fm)
}
show(y | trunc(lb = 1, ub = 30) ~ x, poisson())
show(y | trunc(lb = 1) ~ x, poisson())
show(y | trunc(ub = 30) ~ x, poisson())
show(y | trunc(lb = 0, ub = 30) ~ x, poisson())
show(y | trunc(lb = 1, ub = 30) ~ x, zero_inflated_poisson())
cat("==== the zero_inflated_poisson_lccdf / lcdf functions\n")
sc <- stancode(y | cens(cc, y2) ~ x, data = d, family = zero_inflated_poisson())
cat(sc)
cat("\n==== R-side truncation checks (make_standata)\n")
tr <- function(lab, expr) {
  r <- tryCatch({ force(expr); "ok" }, error = function(e) conditionMessage(e))
  cat(sprintf("%-50s %s\n", lab, r))
}
d0 <- data.frame(y = c(0, 1, 2, 5, 9), x = rnorm(5))
tr("poisson trunc(lb = 0), y has 0", standata(y | trunc(lb = 0) ~ x, d0, poisson()))
tr("poisson trunc(lb = 1), y has 0", standata(y | trunc(lb = 1) ~ x, d0, poisson()))
tr("poisson trunc(ub = 9), y max 9", standata(y | trunc(ub = 9) ~ x, d0, poisson()))
tr("poisson trunc(ub = 8), y max 9", standata(y | trunc(ub = 8) ~ x, d0, poisson()))
tr("poisson trunc(lb = 2, ub = 2)", standata(y | trunc(lb = 2, ub = 2) ~ x, d0, poisson()))
d0$lbv <- c(-Inf, 0, 2, 3, 9); d0$ubv <- c(Inf, 5, 2, 5, 9)
tr("poisson per-row lb/ub incl lb == ub == y",
   standata(y | trunc(lb = lbv, ub = ubv) ~ x, d0, poisson()))
tr("gaussian trunc(lb = -Inf, ub = Inf)",
   standata(y | trunc(lb = -Inf, ub = Inf) ~ x, d0, gaussian()))
d0$lbn <- c(NA, 0, 0, 0, 0)
tr("gaussian trunc(lb = NA row)", standata(y | trunc(lb = lbn) ~ x, d0, gaussian()))
sd0 <- standata(y | trunc(lb = lbv, ub = ubv) ~ x, d0, poisson())
cat("lb:", sd0$lb, " ub:", sd0$ub, "\n")
sd1 <- standata(y | trunc(ub = 9) ~ x, d0, poisson())
cat("poisson trunc(ub = 9) standata lb:", sd1$lb, "\n")
sc <- stancode(y | trunc(ub = 9) ~ x, d0, poisson())
cat(grep("lcdf|lccdf|T\\[", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
