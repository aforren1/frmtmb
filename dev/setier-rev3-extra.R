# Reviewer of lane setier, final check: dev/setier-rev3-short.R part 1
# fits that get the boundary message although lme4 does not call them
# singular. logLik against lme4 and lme4's variance components.
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
src <- readLines("dev/setier-rev3-short.R")
eval(parse(text = src[grep("^gen <- function", src):
                        (grep("^lme <- function", src) + 3L)]))
for (des in c("bin60", "both0")) for (s in 1:40) {
  g <- gen(des, s)
  m <- character()
  f <- withCallingHandlers(frm(g$ff, family = g$fam, data = g$d),
    warning = function(x) invokeRestart("muffleWarning"),
    message = function(x) {m <<- c(m, conditionMessage(x))
      invokeRestart("muffleMessage")})
  l4 <- lme(g)
  if (!any(grepl("^Boundary", m)) || isSingular(l4)) next
  vc <- as.data.frame(VarCorr(l4))
  cat(sprintf("%s seed %d: logLik - lme4 %.3g | lme4 %s | message: %s\n",
              des, s, as.numeric(logLik(f)) - as.numeric(logLik(l4)),
              paste(signif(vc$sdcor, 3), collapse = " "),
              substr(m[grepl("^Boundary", m)][1], 1, 90)))
}
