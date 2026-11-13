# Lane names, backlog 3164: frmtmb on the brms cases of
# dev/names-brms-3164.R. Usage: Rscript dev/names-frm-3164.R <base|lane>
# Output: dev/names-out/frm-3164-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-names-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
sink(sprintf("dev/names-out/frm-3164-%s.txt", arm), split = TRUE)
cat("lib:", find.package("frmtmb"), "\n")
set.seed(3164)
n <- 90
d <- data.frame(x = rnorm(n))
d$yc <- factor(sample(c("a", "b", "c"), n, TRUE),
               levels = c("a", "b", "c", "unused"))
show <- function(lab, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(lab, ":", if (is.character(r)) r else variables(r), "\n")
  invisible(r)
}
show("categorical, unused level", frm(bf(yc ~ x), family = categorical(),
                                       data = d))
show("categorical(K = 3), factor a b c",
     frm(bf(yc ~ x), family = categorical(K = 3), data = droplevels(d)))
d$tr <- 10
d$y1 <- rbinom(n, 10, 0.3)
d$y2 <- rbinom(n, 10 - d$y1, 0.5)
d$y3 <- d$tr - d$y1 - d$y2
show("multinomial(), named cbind",
     frm(bf(cbind(y1, y2, y3) | trials(tr) ~ x), family = multinomial(),
         data = d))
show("multinomial(K = 3), named cbind",
     frm(bf(cbind(y1, y2, y3) | trials(tr) ~ x), family = multinomial(K = 3),
         data = d))
fit <- suppressMessages(frm(bf(yc ~ x), family = categorical(),
                            data = droplevels(d)))
ds <- suppressWarnings(suppressMessages(frm_sample(
  fit, chains = 1, iter = 60, warmup = 30, seed = 1, refresh = 0)))
pe <- posterior_epred(ds)
cat("draws posterior_epred dims:", dim(pe), "\n")
str(dimnames(pe))
cat("draws posterior_predict dimnames:\n")
str(dimnames(posterior_predict(ds)))
sink()
