# Reviewer of lane optima, re-check m5: frmtmb_control(mo_search = FALSE)
# turns the search off and leaves the chart; the 200-seed study's
# seeds 1 to 40 with and without it.
#   Rscript dev/optima-rev2-mosearch.R
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
dat <- function(seed) {
  set.seed(seed)
  lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
  income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
  ls <- c(30, 60, 70, 75)[income] + rnorm(100, sd = 7)
  d <- data.frame(income, ls)
  d$age <- rnorm(100, mean = 40, sd = 10)
  d
}
ev <- c(on = 0, off = 0)
ran <- c(on = 0, off = 0)
worse <- 0
for (s in 1:40) {
  a <- suppressWarnings(frm(ls ~ mo(income) * age, data = dat(s)))
  b <- suppressWarnings(frm(ls ~ mo(income) * age, data = dat(s),
                            control = frmtmb_control(mo_search = FALSE)))
  ev <- ev + c(a$opt$evals, b$opt$evals)
  ran <- ran + c(!is.null(a$opt$mo_search), !is.null(b$opt$mo_search))
  if (as.numeric(logLik(b)) < as.numeric(logLik(a)) - 1e-6) worse <- worse + 1
}
cat("search records: on", ran[["on"]], "off", ran[["off"]],
    "| evaluations: on", ev[["on"]], "off", ev[["off"]],
    "| off below on by > 1e-6:", worse, "of 40\n")
r <- tryCatch(frmtmb_control(mo_search = "no"), error = function(e) e)
cat("mo_search = \"no\":", if (inherits(r, "error")) conditionMessage(r) else
  "accepted", "\n")
