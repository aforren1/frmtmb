# Reviewer, claim 4: read brms 2.23.0's hurdle_cumulative Stan functions
# and R-side log_lik directly (no compile), and evaluate brms's R-side
# log_lik against frmtmb's density on a disc-modeled fit.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(18)
n <- 300
d <- data.frame(x = rnorm(n))
u <- stats::rlogis(n) + 0.8 * d$x
yc <- 1L + (u > -1) + (u > 0.3) + (u > 1.5)
d$y <- ifelse(runif(n) < plogis(-0.6 + 0.5 * d$x), 0L, yc)
show_fun <- function(code, name) {
  lines <- strsplit(code, "\n")[[1]]
  st <- grep(paste0("real ", name, "\\("), lines)
  if (!length(st)) { cat("  (", name, "not in the program)\n"); return() }
  en <- st[1] + 40
  cat(lines[st[1]:min(en, length(lines))], sep = "\n")
}
for (link in c("logit", "probit")) {
  code <- stancode(bf(y ~ x, disc ~ 0 + x), data = d,
                   family = hurdle_cumulative(link))
  cat("\n==== link", link, "disc ~ 0 + x ====\n")
  cat("target line:", grep("target \\+=", strsplit(code, "\n")[[1]],
                           value = TRUE)[1], "\n")
  show_fun(code, paste0("hurdle_cumulative_", link, "_lpmf"))
  sd <- standata(bf(y ~ x, disc ~ 0 + x), data = d,
                 family = hurdle_cumulative(link))
  cat("standata nthres", sd$nthres, " range(Y)", range(sd$Y), "\n")
}
cat("\n==== link logit, disc not modeled (builtin path) ====\n")
code <- stancode(y ~ x, data = d, family = hurdle_cumulative("logit"))
cat(grep("target \\+=|ordered_logistic", strsplit(code, "\n")[[1]],
         value = TRUE), sep = "\n")
show_fun(code, "hurdle_cumulative_ordered_logistic_lpmf")
cat("\n==== brms:::stan_hurdle_ordinal_lpmf source ====\n")
print(brms:::stan_hurdle_ordinal_lpmf)
cat("\n==== brms:::log_lik_hurdle_cumulative ====\n")
print(brms:::log_lik_hurdle_cumulative)
cat("\n==== brms:::posterior_epred_hurdle_cumulative ====\n")
print(brms:::posterior_epred_hurdle_cumulative)
