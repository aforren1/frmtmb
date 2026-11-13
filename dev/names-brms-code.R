# Lane names: brms 2.23.0's Stan declaration order for the predictor
# blocks and group ids, which reorder_pars() keeps within a class
# (no compile). Output: dev/names-out/brms-code.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
source("dev/names-cases.R")
cases <- names_cases(brms::bf, list(
  categorical = categorical(), skew_normal = skew_normal(),
  gaussian = gaussian(), poisson = poisson()))
sink("dev/names-out/brms-code.txt", split = TRUE)
show <- function(nm, f, fam, d) {
  cat("\n=====", nm, "=====\n")
  sc <- strsplit(stancode(f, data = d, family = fam), "\n")[[1]]
  i <- grep("^(parameters|transformed parameters|generated quantities)",
            sc)
  keep <- grepl("^\\s*(vector|real|array|matrix|simplex|cholesky|corr)", sc)
  sec <- cumsum(seq_along(sc) %in% i)
  cat(paste(sec[keep & sec > 0], trimws(sc[keep & sec > 0])), sep = "\n")
  re <- brms:::tidy_ranef(brmsterms(f), data = d)
  print(re[, c("id", "group", "gn", "coef", "cn", "nlpar", "dpar",
               if ("resp" %in% names(re)) "resp")])
}
for (nm in c("skew", "dist", "mv", "me", "gpby", "latres", "twoblk")) {
  cs <- cases[[nm]]
  tryCatch(show(nm, cs$f, cs$fam, cs$d), error = function(e)
    cat("ERROR", conditionMessage(e), "\n"))
}
cat("\n===== skew, formula written sigma before alpha =====\n")
d <- cases$skew$d
tryCatch(show("skew2", bf(ys ~ (1 | h), sigma ~ (1 | ga), alpha ~ (1 | g)),
              skew_normal(), d), error = function(e)
  cat("ERROR", conditionMessage(e), "\n"))
sink()
