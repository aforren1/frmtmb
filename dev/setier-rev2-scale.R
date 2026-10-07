# Reviewer of lane setier, re-check: is the boundary verdict (se_edge_sd's
# curvature screen of 0.1 and se_at_edge's 1e-6) invariant to the units
# of the response? dev/setier-singular.R's ri20 design and
# dev/setier-rev-smallsd.R's s1 design with y multiplied by 1e-3, 1, 1e3.
#   Rscript dev/setier-rev2-scale.R <lib> [seeds]
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
S <- if (length(args) > 1) as.integer(args[2]) else 40L
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", find.package("frmtmb"), "\n")
gen <- function(des, s) {
  set.seed(s)
  if (des == "ri20") {
    d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
    d$y <- 1 + 0.5 * d$x + rnorm(100)
  } else {
    d <- data.frame(g = factor(rep(1:100, each = 10)), x = rnorm(1000))
    d$y <- 1 + 0.5 * d$x + rnorm(100, 0, 0.15)[d$g] + rnorm(1000)
  }
  d
}
for (des in c("ri20", "s1")) {
  flags <- list()
  for (sc in c(1e-3, 1, 1e3)) {
    fl <- logical(S); sing <- logical(S)
    for (s in seq_len(S)) {
      d <- gen(des, s)
      d$y <- d$y * sc
      m <- character()
      f <- withCallingHandlers(frm(y ~ x + (1 | g), data = d),
        warning = function(x) invokeRestart("muffleWarning"),
        message = function(x) {m <<- c(m, conditionMessage(x))
          invokeRestart("muffleMessage")})
      fl[s] <- any(grepl("^Boundary", m))
      sing[s] <- isSingular(suppressMessages(lmer(y ~ x + (1 | g), data = d,
                                                  REML = FALSE)))
    }
    flags[[as.character(sc)]] <- fl
    cat(sprintf("%-5s y x %-6g | lme4 singular %d | flagged %d (on lme4-singular %d, others %d)\n",
                des, sc, sum(sing), sum(fl), sum(fl & sing), sum(fl & !sing)))
  }
  cat(sprintf("%-5s seeds whose verdict differs between scales: %d\n", des,
              sum(!(flags[[1]] == flags[[2]] & flags[[2]] == flags[[3]]))))
}
