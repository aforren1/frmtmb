# Review of lane meme: the same round trip for cor and sd group rows,
# on the lane library and on rellib-r15.
#   Rscript dev/review-meme-grpinh2.R <lane|base>
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r15",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-meme-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-meme/dev/review-meme-shapes.R")
d <- rv_data()
f <- bf(y ~ x + (1 + x | g))
fit0 <- frm(f, data = d, dry_run = "objective")
ts <- default_prior(f, data = d, route = "sample")
print(as.data.frame(ts)[ts$class %in% c("sd", "cor"),
                        c("prior", "class", "coef", "group")])
t1 <- ts
t1$prior[t1$class == "cor" & !nzchar(t1$group)] <- "lkj(4)"
t1$prior[t1$class == "sd" & !nzchar(t1$group) & !nzchar(t1$coef)] <-
  "exponential(3)"
rs <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit0, t1))
cat("class rows edited (cor lkj(4), sd exponential(3)); sampled:",
    vapply(rs$ri$entries, function(e) {
      paste0(e$comp, "[", paste(e$idx, collapse = ","), "]:", e$dist$kind,
             e$dist$eta %||% e$dist$rate %||% "")
    }, ""), "\n")
