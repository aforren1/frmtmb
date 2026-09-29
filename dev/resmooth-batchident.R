# Lane wt-resmooth, nits round. The batching change must not move a
# single number: a whole-block batch perturbs every column of the block
# at once, and a row's value moves only through its own column because
# every other column's design entry on that row is exactly zero, so each
# Jacobian entry is the one a singleton batch would have produced. That
# is an exactness ARGUMENT; this script is the check, by identical() on
# the saved Est.Error arrays of the reviewer's own cells
# (dev/resmooth-rev-estse.R constructions A, B, C, D, F, I, plus the
# gp() and hsgp() cells the fix is for and a no-smooth control).
#   RESMOOTH_ARM=before Rscript dev/resmooth-batchident.R
#   RESMOOTH_ARM=after  Rscript dev/resmooth-batchident.R
#   RESMOOTH_ARM=cmp    Rscript dev/resmooth-batchident.R
arm <- Sys.getenv("RESMOOTH_ARM")
out <- file.path("dev", paste0("resmooth-batchident-", arm, ".rds"))
if (identical(arm, "cmp")) {
  a <- readRDS("dev/resmooth-batchident-before.rds")
  b <- readRDS("dev/resmooth-batchident-after.rds")
  stopifnot(identical(names(a), names(b)))
  ok <- TRUE
  for (nm in names(a)) {
    id <- identical(a[[nm]], b[[nm]])
    mx <- if (is.numeric(a[[nm]]) && is.numeric(b[[nm]])) {
      max(abs(as.vector(a[[nm]]) - as.vector(b[[nm]])))
    } else NA_real_
    cat(sprintf("%-46s identical %-5s max abs diff %s\n", nm, id,
                if (is.na(mx)) "n/a" else sprintf("%.3e", mx)))
    ok <- ok && id
  }
  cat(if (ok) "ALL IDENTICAL\n" else "A CELL MOVED\n")
  quit(status = if (ok) 0L else 1L)
}
lib <- switch(arm,
              before = "C:/Users/adf44/source/r/wt-resmooth-lib",
              after = "C:/Users/adf44/source/r/wt-resmooth-lib2",
              stop("set RESMOOTH_ARM to before, after or cmp"))
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("ARM:", arm, "from", find.package("frmtmb"), "\n")
res <- list()
take <- function(lab, fit, newdata) {
  for (rf in list(NA, NULL)) {
    key <- paste0(lab, " re_formula=", if (is.null(rf)) "NULL" else "NA")
    v <- tryCatch(fitted(fit, newdata = newdata, re_formula = rf),
                  error = function(e) paste("ERROR:",
                                            conditionMessage(e)))
    res[[key]] <<- v
    k2 <- paste0(lab, " in sample re_formula=",
                 if (is.null(rf)) "NULL" else "NA")
    res[[k2]] <<- tryCatch(fitted(fit, re_formula = rf),
                           error = function(e) paste("ERROR:",
                                                     conditionMessage(e)))
  }
}

set.seed(23)
n <- 240
dA <- data.frame(x = stats::runif(n, -2, 2))
latA <- 1.2 * sin(2 * dA$x) + stats::rlogis(n)
dA$y <- factor(cut(latA, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fA <- suppressWarnings(frm(bf(y ~ s(x, k = 8)), family = cumulative(),
                           data = dA))
take("A cumulative s(x,k=8)", fA, data.frame(x = c(-1.5, 0, 1.5)))

set.seed(31)
ng <- 8
per <- 30
dB <- data.frame(g = factor(rep(seq_len(ng), each = per)),
                 x = stats::runif(ng * per, -2, 2))
latB <- 1.0 * sin(1.5 * dB$x) + stats::rnorm(ng, 0, 0.8)[dB$g] +
  stats::rlogis(nrow(dB))
dB$y <- factor(cut(latB, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fB <- suppressWarnings(frm(bf(y ~ s(x, k = 8) + (1 | g)),
                           family = cumulative(), data = dB))
take("B cumulative s(x,k=8)+(1|g)", fB, dB[c(5, 40, 200), c("x", "g")])

set.seed(37)
nC <- 250
dC <- data.frame(x = stats::runif(nC, -2, 2), z = stats::runif(nC, -2, 2))
latC <- sin(dC$x) + 0.6 * dC$z + stats::rlogis(nC)
dC$y <- factor(cut(latC, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fC <- suppressWarnings(frm(bf(y ~ t2(x, z, k = 4)), family = cumulative(),
                           data = dC))
take("C cumulative t2(x,z,k=4)", fC,
     data.frame(x = c(-1, 0, 1), z = c(0.5, -0.5, 1)))

set.seed(41)
nD <- 160
dD <- data.frame(x = stats::runif(nD, -2, 2))
latD <- 1.1 * sin(2 * dD$x) + stats::rlogis(nD)
dD$y <- factor(cut(latD, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fD <- suppressWarnings(frm(bf(y ~ gp(x)), family = cumulative(), data = dD))
take("D cumulative gp(x) n_b=160", fD, data.frame(x = c(-1.5, 0, 1.5)))
# the cell the fix is for: three OBSERVED positions, indicator design
take("D2 cumulative gp(x) observed rows", fD, dD[c(1, 5, 9), "x",
                                                 drop = FALSE])
fG <- suppressWarnings(frm(bf(y ~ gp(x, k = 12, c = 5 / 4)),
                          family = cumulative(), data = dD))
take("D3 cumulative hsgp(x) k=12", fG, data.frame(x = c(-1.5, 0, 1.5)))

set.seed(43)
ngF <- 4
perF <- 70
dF <- data.frame(g = factor(rep(seq_len(ngF), each = perF)),
                 x = stats::runif(ngF * perF, -2, 2))
latF <- sin(1.5 * dF$x) + stats::rlogis(nrow(dF))
cutsF <- list(c(-1.5, 0.5), c(-0.5, 1.0), c(-1.0, 0.0), c(-0.2, 1.4))
dF$y <- ordered(vapply(seq_len(nrow(dF)), function(i) {
  cc <- cutsF[[as.integer(dF$g[i])]]
  1L + sum(latF[i] > cc)
}, 1L))
fF <- suppressWarnings(
  frm(bf(y | thres(gr = g) ~ s(x, k = 6)), family = cumulative(), data = dF))
take("F cumulative thres(gr=g) s(x,k=6)", fF,
     dF[c(3, 80, 150, 250), c("x", "g")])

set.seed(47)
nI <- 400
dI <- data.frame(x = stats::runif(nI, -2, 2))
latI <- 0.9 * dI$x + stats::rlogis(nI)
dI$y <- factor(cut(latI, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fI <- suppressWarnings(frm(bf(y ~ s(x, k = 8)), family = cumulative(),
                           data = dI))
take("I cumulative s(x,k=8) near-linear", fI,
     data.frame(x = c(-1.5, 0, 1.5)))

# a (1 + x | g) block, whose whole-block test MUST fail so that the
# per-position split still runs: the case the change could break
set.seed(67)
ngJ <- 12
perJ <- 25
dJ <- data.frame(g = factor(rep(seq_len(ngJ), each = perJ)),
                 x = stats::runif(ngJ * perJ, -2, 2))
bJ <- matrix(stats::rnorm(ngJ * 2, 0, 0.6), ngJ, 2)
latJ <- sin(dJ$x) + bJ[dJ$g, 1L] + bJ[dJ$g, 2L] * dJ$x +
  stats::rlogis(nrow(dJ))
dJ$y <- factor(cut(latJ, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fJ <- suppressWarnings(frm(bf(y ~ s(x, k = 6) + (1 + x | g)),
                           family = cumulative(), data = dJ))
take("J cumulative s(x,k=6)+(1+x|g)", fJ, dJ[c(2, 30, 90), c("x", "g")])

# a multi-membership block, whose batching returns NULL on both arms
set.seed(71)
nM <- 200
dM <- data.frame(x = stats::runif(nM, -2, 2),
                 g1 = factor(sample(letters[1:6], nM, TRUE)),
                 g2 = factor(sample(letters[1:6], nM, TRUE)))
latM <- sin(dM$x) + stats::rlogis(nM)
dM$y <- factor(cut(latM, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fM <- tryCatch(suppressWarnings(
  frm(bf(y ~ s(x, k = 6) + (1 | mm(g1, g2))), family = cumulative(),
      data = dM)), error = function(e) NULL)
if (!is.null(fM)) {
  take("M cumulative s(x,k=6)+(1|mm(g1,g2))", fM,
       dM[c(1, 50, 120), c("x", "g1", "g2")])
}

# CONTROL: no smooth, no gp, so nothing in this round can reach it
set.seed(59)
dH <- data.frame(g = factor(rep(seq_len(40), each = 20)),
                 x = stats::runif(800, -2, 2))
latH <- 0.8 * dH$x + stats::rnorm(40, 0, 0.7)[dH$g] +
  stats::rlogis(nrow(dH))
dH$y <- factor(cut(latH, c(-Inf, -0.8, 0.8, Inf), labels = FALSE),
               ordered = TRUE)
fH <- suppressWarnings(frm(bf(y ~ x + (1 | g)), family = cumulative(),
                           data = dH))
take("CONTROL cumulative x+(1|g)", fH, dH[c(1, 300, 700), c("x", "g")])

saveRDS(res, out)
cat("saved", out, "with", length(res), "cells\n")
