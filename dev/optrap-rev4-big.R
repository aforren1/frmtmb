# Reviewer of lane optrap, consolidation check: the designs that turn on
# newton_polish() for a fit with random effects by themselves
# (dev/optrap-rev4-re.R's second part), on any build, against lmer.
# year_sq: y ~ year + I(year^2) + (1 | g), year 1990 to 2020; big:
# y ~ big + (1 | g), big = 1e5 + U(0, 1). Seeds 1 to 10.
#   Rscript dev/optrap-rev4-big.R <lib or none> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
lab <- a[2]
cat(lab, "lib", dirname(find.package("frmtmb")), "\n")
for (s in 1:10) {
  set.seed(s)
  n <- 400
  d <- data.frame(g = factor(rep(1:40, each = 10)),
                  year = sample(1990:2020, n, TRUE), big = 1e5 + runif(n))
  d$y <- 0.002 * (d$year - 2005)^2 + 0.1 * (d$year - 2005) +
    rnorm(40, 0, 0.5)[d$g] + rnorm(n)
  d$y2 <- 0.5 * (d$big - 1e5) + rnorm(40, 0, 0.5)[d$g] + rnorm(n)
  for (des in list(list("year_sq", y ~ year + I(year^2) + (1 | g)),
                   list("big", y2 ~ big + (1 | g)))) {
    w <- character()
    f <- withCallingHandlers(suppressMessages(frm(des[[2]], data = d)),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    l4 <- suppressMessages(suppressWarnings(lmer(des[[2]], data = d,
                                                 REML = FALSE)))
    cat(sprintf("%s seed %2d %-7s | code %d | logLik - lmer %9.3g | steps %s | %s\n",
                lab, s, des[[1]], f$opt$convergence,
                as.numeric(logLik(f)) - as.numeric(logLik(l4)),
                f$opt$polish[["steps"]] %||% 0,
                if (length(w)) substr(paste(substr(w, 1, 45), collapse = " | "),
                                      1, 120) else "silent"))
  }
}
