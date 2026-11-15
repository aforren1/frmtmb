# What a fit keeps of the rows it read: fit$data, frame$y, subset_rows.
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(6)
du <- data.frame(x = rnorm(12), z = rnorm(12))
du$y <- du$x + rnorm(12)
du$x[1:2] <- NA
ux <- suppressMessages(frm(y ~ x, data = du))
str(ux$data)
cat("rownames(fit$data):", rownames(ux$data), "\n")
cat("n_obs", ux$frame$n_obs, " y len", length(ux$frame$y$y), "\n")
print(ux$frame$na_action)
d6 <- data.frame(x = rnorm(10), y1 = rnorm(10), y2 = rnorm(10),
                 s = rep(c(TRUE, FALSE), 5))
f <- frm(bf(y1 | subset(s) ~ x) + bf(y2 ~ x) + set_rescor(FALSE), data = d6)
str(f$frame$y)
str(f$frame$subset_rows)
cat("rownames(fit$data):", rownames(f$data), "\n")
names(f$frame)
