# Reviewer, lane wt-conditions (recheck 1): tabulate
# dev/conditions-rev-newdata-brms.R, base against lane against brms.
#   Rscript dev/conditions-rev-newdata-brms-cmp.R
o <- "C:/Users/adf44/source/r/frmtmb-wt-conditions/dev/conditions-rev-log"
b <- readRDS(file.path(o, "newdata-brms-base.rds"))
l <- readRDS(file.path(o, "newdata-brms-lane.rds"))
stopifnot(identical(b[, 1:3], l[, 1:3]), identical(b$brms, l$brms))
st <- function(x) ifelse(x == "OK", "OK", "ERR")
tab <- data.frame(b[, 1:3], base = st(b$frmtmb), lane = st(l$frmtmb),
                  brms = st(b$brms))
cat("rows", nrow(tab), "\n")
cat("base OK, lane ERR (regression):", sum(tab$base == "OK" & tab$lane == "ERR"), "\n")
cat("base ERR, lane OK:", sum(tab$base == "ERR" & tab$lane == "OK"), "\n")
cat("lane ERR where brms OK:", sum(tab$lane == "ERR" & tab$brms == "OK"), "\n")
cat("lane OK where brms ERR:", sum(tab$lane == "OK" & tab$brms == "ERR"), "\n\n")
options(width = 250)
dis <- tab$base != tab$lane | tab$lane != tab$brms
x <- cbind(tab, lane_msg = substr(l$frmtmb, 1, 75), brms_msg = substr(b$brms, 1, 60))
print(x[dis, ], row.names = FALSE, right = FALSE)
