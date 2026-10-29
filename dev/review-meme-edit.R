# Review of lane meme: a sample-route default table, its class row
# edited, passed back. The coef rows print "(flat)"; do they inherit the
# edited class row (brms's "" does) or override it?
#   Rscript dev/review-meme-edit.R > dev/review-meme-out/edit.txt
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source("C:/Users/adf44/source/r/frmtmb-wt-meme/dev/review-meme-shapes.R")
d <- rv_data()
f <- bf(y ~ me(x, sx) + me(z, sz) + w + (1 + w | g))
fit0 <- frm(f, data = d, dry_run = "objective")
show <- function(lab, tab) {
  rs <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit0, tab))
  cat("--", lab, "\n  ")
  cat(vapply(rs$ri$entries, function(e) {
    paste0(e$comp, "[", paste(e$idx, collapse = ","), "]:", e$dist$kind,
           if (!is.null(e$dist$rate)) e$dist$rate,
           if (!is.null(e$dist$eta)) e$dist$eta)
  }, ""), fill = 76)
}
ts <- default_prior(f, data = d, route = "sample")
show("sample-route table as it is", ts)
t1 <- ts
t1$prior[t1$class == "sdme" & !nzchar(t1$coef)] <- "exponential(3)"
show("sdme class row edited to exponential(3)", t1)
t2 <- ts
t2$prior[t2$class == "sd" & !nzchar(t2$coef) & !nzchar(t2$group)] <-
  "exponential(3)"
t2$prior[t2$class == "sd" & !nzchar(t2$coef) & nzchar(t2$group)] <-
  "exponential(3)"
show("sd class and group rows edited to exponential(3)", t2)
t3 <- ts
t3$prior[t3$class == "meanme" & !nzchar(t3$coef)] <- "normal(0, 3)"
show("meanme class row edited to normal(0, 3)", t3)
t4 <- ts
t4$prior[t4$class == "corme"] <- "lkj(4)"
show("corme row edited to lkj(4)", t4)
cat("DONE\n")
