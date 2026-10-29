# Review of lane meme: a negative or positive lb on sdme, and on sd, on
# both routes. Rscript dev/review-meme-negsd.R > dev/review-meme-out/negsd.txt
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
q <- function(lab, expr) {
  r <- tryCatch(suppressMessages(suppressWarnings(expr)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("--", lab, "\n")
  print(r)
}
f <- bf(y ~ me(x, sx) + me(z, sz))
for (lb in c(-1, 0.2)) {
  pr <- set_prior("exponential(1)", class = "sdme", lb = lb)
  q(paste("frm sdme lb", lb), {
    ft <- frm(f, data = d, prior = pr)
    ri <- frmtmb:::resolve_priorlist(ft, pr)
    list(lower = ri$lower[grep("logsdme", names(ri$lower))],
         sdme = exp(ft$estimates$logsdme))
  })
  q(paste("frm_sample sdme lb", lb), {
    ds <- frm_sample(f, data = d, family = gaussian(), prior = pr,
                     chains = 1, iter = 300, seed = 2, refresh = 0)
    range(as.matrix(ds)[, c("sdme_mex", "sdme_mez")])
  })
}
fr <- bf(y ~ x + (1 | g))
q("frm sd lb -1", {
  pr <- set_prior("exponential(1)", class = "sd", lb = -1)
  ft <- frm(fr, data = d, prior = pr)
  frmtmb:::resolve_priorlist(ft, pr)$lower
})
cat("DONE\n")
