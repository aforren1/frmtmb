# An ORDERED factor inside cs(): model.matrix() gives contr.poly
# columns, so does frmtmb match brms there too?
#   Rscript dev/csfactor-ordered.R > dev/csfactor-log/ordered.txt
.libPaths(c("C:/Users/adf44/source/r/wt-csfactor-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(405)
n <- 400
fo <- factor(sample(c("lo", "mid", "hi"), n, TRUE),
             levels = c("lo", "mid", "hi"), ordered = TRUE)
eff <- c(lo = -1, mid = 0, hi = 1.5)[as.character(fo)]
p1 <- plogis(-0.3 + eff); p2 <- (1 - p1) * plogis(0.5 - eff)
u <- runif(n)
d <- data.frame(fo = fo,
                yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
ff <- frm(bf(yo ~ cs(fo)), family = sratio(), data = d)
cat("frmtmb variables():",
    paste(grep("^bcs", variables(ff), value = TRUE), collapse = " "), "\n")
cat("frmtmb logLik:", sprintf("%.6f", as.numeric(logLik(ff))), "\n")
cat("frmtmb fitted at the three levels:\n")
print(round(frm_linpred(ff, newdata = data.frame(
  fo = factor(c("lo", "mid", "hi"),
              levels = c("lo", "mid", "hi"), ordered = TRUE)),
  type = "response"), 6))
cat("empirical shares:\n")
print(round(prop.table(table(d$fo, d$yo), 1), 6))
# a one-row newdata keeping only the level it needs
cat("one row, fo = ordered('hi') with all three levels declared:\n")
print(round(frm_linpred(ff, newdata = data.frame(
  fo = factor("hi", levels = c("lo", "mid", "hi"), ordered = TRUE)),
  type = "response"), 6))
sd <- brms::standata(brms::bf(yo ~ cs(fo)), family = brms::sratio(),
                     data = d)
cat("brms Kcs =", sd$Kcs, "Xcs columns:",
    paste(colnames(sd$Xcs), collapse = ", "), "\n")
cat("brms Xcs head:\n")
print(utils::head(sd$Xcs, 3))
cat("frmtmb stored cs columns:\n")
fr <- frm(bf(yo ~ cs(fo)), family = sratio(), data = d,
          dry_run = "frame")
cs <- fr$linpreds[["yo.mu"]][["cs"]]
for (ct in cs) {
  cat(" ", ct[["label"]], ":", paste(round(utils::head(ct[["vals"]], 3), 6),
                                     collapse = " "), "\n")
}
