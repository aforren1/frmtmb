A <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/bitwise-p1base1.rds")
B <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/bitwise-p1lane.rds")
k <- grep("[|]slope[|]", names(A), value = TRUE)
els <- c("error", "ll", "par", "conv", "se_full", "se_fix", "sd_re")
cat("slope fits", length(k), "; errored", sum(vapply(A[k], function(r) !is.null(r$error), NA)),
    "; identical", sum(vapply(k, function(i) identical(A[[i]][els], B[[i]][els]), NA)), "\n")
