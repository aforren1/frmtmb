.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
try_show <- function(label, expr) {
  cat("\n==========", label, "\n")
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)),
                warning = function(w) paste("WARNING:", conditionMessage(w)))
  print(r)
  invisible(r)
}
# --- rate ---------------------------------------------------------------
dr <- data.frame(y = rpois(10, 1), x = rnorm(10), time = 1:10,
                 cc = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0))
cat("\n=== rate poisson\n")
try_show("sc", cat(stancode(y | rate(time) ~ x, dr, poisson())))
cat("\n=== rate poisson identity\n")
try_show("sc", cat(stancode(y | rate(time) ~ x, dr, poisson("identity"))))
cat("\n=== rate negbinomial\n")
try_show("sc", cat(stancode(y | rate(time) ~ x, dr, negbinomial())))
cat("\n=== rate negbinomial shape ~ x\n")
try_show("sc", cat(stancode(bf(y | rate(time) ~ x, shape ~ x), dr, negbinomial())))
cat("\n=== rate negbinomial2\n")
try_show("sc", cat(stancode(y | rate(time) ~ x, dr, brmsfamily("negbinomial2"))))
cat("\n=== rate geometric\n")
try_show("sc", cat(stancode(y | rate(time) ~ x, dr, geometric())))
cat("\n=== rate poisson cens\n")
try_show("sc", cat(stancode(y | rate(time) + cens(cc) ~ x, dr, poisson())))
try_show("rate gaussian", stancode(y | rate(time) ~ x, dr, gaussian()))
try_show("rate zip", stancode(y | rate(time) ~ x, dr, zero_inflated_poisson()))
try_show("rate 0", standata(y | rate(time - 1) ~ x, dr, poisson()))
try_show("rate const", standata(y | rate(2) ~ x, dr, poisson())$denom)

# --- cat --------------------------------------------------------------
dc <- data.frame(s = sample(1:5, 9, TRUE), x = rnorm(9))
try_show("cat standata", {
  w <- NULL
  r <- withCallingHandlers(standata(s | cat(6) ~ x, dc, cumulative()),
                           warning = function(ww) {
                             w <<- c(w, conditionMessage(ww))
                             invokeRestart("muffleWarning")
                           })
  list(nthres = r$nthres, warnings = w)
})
try_show("cat on gaussian", standata(s | cat(6) ~ x, dc, gaussian()))
try_show("cat + thres", standata(s | cat(6) + thres(5) ~ x, dc, cumulative()))
