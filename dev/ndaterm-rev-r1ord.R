# Reviewer of lane ndaterm, fix round 1: NA draws of an ordinal
# mixture's predict(propagate_error = TRUE), and mixture_mvn() on
# newdata. Usage: Rscript dev/ndaterm-rev-r1ord.R <arm>. Seed 7.
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "mine"
base <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
            file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(7)
do <- data.frame(x = rnorm(400))
ko <- rbinom(400, 1, 0.3)
do$yo <- cut(do$x + ifelse(ko == 1, 2, 0) + rlogis(400),
             c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
f4 <- frm(bf(yo ~ x), family = mixture(cumulative(), cumulative()),
          data = do)
ndo <- data.frame(x = c(-1, 0, 1))
for (nd in list(ndo, NULL)) {
  set.seed(13)
  w <- character()
  p <- withCallingHandlers(
    predict(f4, newdata = nd, summary = FALSE, ndraws = 4000),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning")
    })
  pn <- matrix(as.numeric(p), nrow(p))
  cat(if (is.null(nd)) "in sample" else "newdata",
      ": share of NA draws", mean(is.na(pn)), "; warnings:", length(w),
      if (length(w)) substr(w[1], 1, 100), "\n")
}
s <- summary(predict(f4, newdata = ndo, ndraws = 4000))
print(head(predict(f4, newdata = ndo, ndraws = 4000)))
fc <- frm(bf(yo ~ x), family = cumulative(), data = do)
set.seed(13)
p <- predict(fc, newdata = ndo, summary = FALSE, ndraws = 4000)
cat("plain cumulative: share of NA draws",
    mean(is.na(matrix(as.numeric(p), nrow(p)))), "\n")
cat("\n== mixture_mvn(K = 2, D = 2) on newdata\n")
set.seed(7)
n <- 300
k <- rbinom(n, 1, 0.3)
dm <- data.frame(x = rnorm(n))
dm$Y <- cbind(rnorm(n, ifelse(k == 1, 4, 0) + dm$x),
              rnorm(n, ifelse(k == 1, -3, 0)))
fm <- tryCatch(frm(bf(Y ~ x), family = mixture_mvn(K = 2, D = 2),
                   data = dm), error = function(e) conditionMessage(e))
if (is.character(fm)) cat("fit:", fm, "\n") else {
  nd <- data.frame(x = c(-1, 0, 1))
  f <- tryCatch(fitted(fm, newdata = nd),
                error = function(e) conditionMessage(e))
  cat("fitted(newdata):\n"); print(f)
  set.seed(14)
  p <- tryCatch(predict(fm, newdata = nd, summary = FALSE, ndraws = 4000,
                        propagate_error = FALSE),
                error = function(e) conditionMessage(e))
  if (is.character(p)) cat("predict(newdata):", p, "\n") else {
    cat("predict(newdata) dims:", dim(p), "\n")
    print(apply(p, seq_along(dim(p))[-1], mean))
  }
  s <- tryCatch(simulate(fm, nsim = 2000, newdata = nd, seed = 14),
                error = function(e) conditionMessage(e))
  if (is.character(s)) cat("simulate(newdata):", s, "\n") else {
    cat("simulate(newdata) str:\n"); str(s[1:2])
  }
}
cat("done\n")
