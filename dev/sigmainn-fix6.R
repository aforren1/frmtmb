# Fix round 1: move the .get_sigma() block of R/insight.R above the
# comment that belongs to get_residuals.frmtmb_fit().
f <- "R/insight.R"
x <- readLines(f)
a <- grep("^# insight::get_residuals[(][)] has no method", x)
b <- grep("^# insight::get_sigma[(][)] is not generic", x)
e <- grep("^[.]get_sigma[.]frmtmb_fit <- function", x)
e <- e + which(x[(e + 1):length(x)] == "}")[1]
stopifnot(length(a) == 1, length(b) == 1, b > a, x[e + 1] == "")
blk <- x[b:(e + 1)]
x <- x[-(b:(e + 1))]
x <- append(x, blk, after = a - 1)
writeLines(x, f)
