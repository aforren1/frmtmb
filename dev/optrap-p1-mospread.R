# Lane optrap, punch round 1: predictor_col_spread() on the mo() entries
# lane moparity's frame builds for mo(x) * f (mox, mox:fb, mox:fc: one
# placeholder column each, `codes` the codes of x, `mult` NULL for the
# main term and the 0/1 contrast column of mo_sp_design() otherwise),
# and with a factor value missing (model.matrix(na.pass) gives NA).
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
spread <- frmtmb:::predictor_col_spread
set.seed(1)
n <- 120
codes <- sample(0:4, n, TRUE)
f <- factor(sample(c("a", "b", "c"), n, TRUE))
mm <- model.matrix(~ f)
X <- cbind("(Intercept)" = rep(1, n), mox = 0, "mox:fb" = 0, "mox:fc" = 0)
lp <- list(X = X, mo = list(
  list(col = 2L, codes = codes, mult = NULL),
  list(col = 3L, codes = codes, mult = unname(mm[, "fb"])),
  list(col = 4L, codes = codes, mult = unname(mm[, "fc"]))))
cat("spreads:", spread(lp), "\n")
cat("by hand:", 2^round(log2(c(sd(codes), sd(codes * mm[, "fb"]),
                               sd(codes * mm[, "fc"])))), "\n")
lp$mo[[2]]$mult[5] <- NA
cat("with an NA in mox:fb's multiplier:", spread(lp), "\n")
lp$mo[[3]]$mult <- lp$mo[[3]]$mult[-1]
cat("with mox:fc's multiplier one row short:", spread(lp), "\n")
