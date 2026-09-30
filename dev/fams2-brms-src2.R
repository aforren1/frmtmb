# How brms 2.23.0 treats the hurdle_cumulative response and its extra
# category: grep every function in the namespace for "extra_cat".
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
ns <- asNamespace("brms")
hit <- character()
for (nm in ls(ns, all.names = TRUE)) {
  f <- get(nm, ns)
  if (!is.function(f)) next
  src <- paste(deparse(f), collapse = "\n")
  if (grepl("extra_cat|hurdle_cumulative", src)) hit <- c(hit, nm)
}
print(hit)
for (nm in hit) {
  cat("\n=====", nm, "=====\n")
  print(get(nm, ns))
}
cat("\n===== data_response.brmsterms ... (lines with extra_cat or ordinal) =====\n")
for (nm in c("data_response.brmsterms", "validate_resp_ordinal",
             "check_discrete_trunc_bounds", "get_thres", "is_ordinal",
             "has_extra_cat", "has_thres", "subset_thres", "pordinal",
             "log_cdf", "log_ccdf", "inv_link_cumulative", "has_cat",
             "has_ordered_thres", "has_equidistant_thres",
             "has_sum_to_zero_thres", "stan_thres", ".stan_thres")) {
  cat("\n=====", nm, "=====\n")
  if (exists(nm, ns)) print(get(nm, ns)) else cat("(absent)\n")
}
# response coding with a factor and with integers
set.seed(2)
d <- data.frame(x = rnorm(30))
d$yf <- factor(sample(c("none", "low", "mid", "high"), 30, TRUE),
               levels = c("none", "low", "mid", "high"), ordered = TRUE)
sd <- tryCatch(standata(bf(yf ~ x), data = d, family = hurdle_cumulative()),
               error = function(e) conditionMessage(e))
cat("\nfactor response:\n"); print(if (is.list(sd)) list(Y = sd$Y, nthres = sd$nthres, levels = table(d$yf)) else sd)
d$yi <- sample(1:4, 30, TRUE)
sd <- tryCatch(standata(bf(yi ~ x), data = d, family = hurdle_cumulative()),
               error = function(e) conditionMessage(e))
cat("\nno zeros:\n"); print(if (is.list(sd)) list(Y = sd$Y, nthres = sd$nthres) else sd)
d$yn <- sample(-1:3, 30, TRUE)
sd <- tryCatch(standata(bf(yn ~ x), data = d, family = hurdle_cumulative()),
               error = function(e) conditionMessage(e))
cat("\nnegative:\n"); print(if (is.list(sd)) list(Y = sd$Y) else sd)
d$yt <- sample(0:3, 30, TRUE)
sd <- tryCatch(standata(bf(yt | thres(5) ~ x), data = d, family = hurdle_cumulative()),
               error = function(e) conditionMessage(e))
cat("\nthres(5):\n"); print(if (is.list(sd)) list(Y = sd$Y, nthres = sd$nthres) else sd)
