# Reviewer 2, item 5: lcdf and lccdf now run only on the censored rows
# (ddm_on_rows()). The objective, evaluated WITHOUT optimizing
# (dry_run = "objective"), at three fixed parameter vectors, on designs
# with 0, some (non-contiguous, mixed codes), and 100 percent censored
# rows. Run once per library; r2-rows-compare.R compares.
# Usage: Rscript r2-rows.R <tag> <lib>
source("dev/phase3b-review2/r2-prelude.R")
a <- commandArgs(trailingOnly = TRUE)
tag <- a[[1]]
.libPaths(c(a[[2]], "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
cat("frmtmb.eam from", find.package("frmtmb.eam"), "\n")

mk <- function(seed, n = 400) {
  set.seed(seed)
  d <- ddm_simulate(n, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
  d$g <- factor(sample(1:8, n, TRUE))
  d$code <- 0L; d$y2 <- d$rt
  d
}
designs <- list()
d <- mk(1); designs$none_cens0 <- list(d = d, f = "cens")
d <- mk(1); designs$none_plain <- list(d = d, f = "plain")
# some: right beyond 1.1 s, left below 0.45 s, interval on a random
# scatter of rows; codes interleave
d <- mk(2)
r <- d$rt > 1.1; d$code[r] <- 1L; d$rt[r] <- 1.1
l <- d$rt < 0.45 & d$code == 0L; d$code[l] <- -1L; d$rt[l] <- 0.45
pick <- which(d$code == 0L)[sample.int(sum(d$code == 0L), 60)]
lo <- floor(d$rt[pick] * 10) / 10
d$code[pick] <- 2L; d$y2[pick] <- lo + 0.1; d$rt[pick] <- pmax(lo, 0.45)
designs$some_mixed <- list(d = d, f = "cens2")
# some, right only, scattered (every third row plus a random set)
d <- mk(3)
k <- sort(unique(c(seq(1, 400, 3), sample.int(400, 30))))
d$code[k] <- 1L; d$rt[k] <- pmax(d$rt[k], 0.9)
designs$some_right <- list(d = d, f = "cens")
# 100 percent right
d <- mk(4); d$code[] <- 1L; d$rt <- pmax(d$rt, 0.8)
designs$all_right <- list(d = d, f = "cens")
# 100 percent left
d <- mk(5); d$code[] <- -1L; d$rt <- pmax(d$rt, 0.6)
designs$all_left <- list(d = d, f = "cens")
# 100 percent interval
d <- mk(6); lo <- floor(d$rt * 10) / 10
d$code[] <- 2L; d$y2 <- lo + 0.1; d$rt <- pmax(lo, 0.35)
designs$all_int <- list(d = d, f = "cens2")
# right censoring under trunc(ub = 3), with a random effect
d <- mk(7); d <- d[d$rt < 3, ]
r <- d$rt > 1.2; d$code[r] <- 1L; d$rt[r] <- 1.2
designs$trunc_right <- list(d = d, f = "trunc")

fit_obj <- function(x) {
  f <- switch(x$f,
    plain = bf(rt | dec(upper) ~ 1 + (1 | g), bs ~ 1, ndt ~ 1, bias = 0.5),
    cens = bf(rt | dec(upper) + cens(code) ~ 1 + (1 | g), bs ~ 1, ndt ~ 1,
              bias = 0.5),
    cens2 = bf(rt | dec(upper) + cens(code, y2) ~ 1 + (1 | g), bs ~ 1,
               ndt ~ 1, bias = 0.5),
    trunc = bf(rt | dec(upper) + cens(code) + trunc(ub = 3) ~ 1 + (1 | g),
               bs ~ 1, ndt ~ 1, bias = 0.5))
  tryCatch(frm(f, family = wiener(), data = x$d, dry_run = "objective"),
           error = function(e) conditionMessage(e))
}
out <- list()
for (nm in names(designs)) {
  o <- fit_obj(designs[[nm]])
  if (is.character(o)) { out[[nm]] <- list(error = o); cat(nm, "ERROR", o, "\n"); next }
  p0 <- o$obj$par
  P <- rbind(p0, p0 + c(0.3, -0.1, 0.4, rep(0.2, length(p0) - 3)),
             p0 + c(-0.5, 0.2, -1, rep(-0.3, length(p0) - 3)))
  fn <- apply(P, 1, function(p) o$obj$fn(p))
  gr <- t(apply(P, 1, function(p) o$obj$gr(p)))
  out[[nm]] <- list(par = P, fn = fn, gr = gr, n = nrow(designs[[nm]]$d),
                    codes = table(designs[[nm]]$d$code))
  cat(sprintf("%-12s n %d  fn %s\n", nm, nrow(designs[[nm]]$d),
              paste(sprintf("%.15g", fn), collapse = "  ")))
}
saveRDS(out, sprintf("dev/phase3b-review2/rows-%s.rds", tag))
