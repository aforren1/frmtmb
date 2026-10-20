# Lane adhess: compare dev/adhess-sfexact.R's two arms row by row.
#   Rscript dev/adhess-sfexact-cmp.R <dir>
d <- commandArgs(TRUE)[1]
rd <- function(arm) do.call(rbind, lapply(c("lin", "poly"), function(g) {
  x <- readRDS(file.path(d, sprintf("%s-%s.rds", g, arm)))
  x$group <- g
  x
}))
L <- rd("lane")
B <- rd("base")
k <- c("group", "model", "seed")
m <- merge(L, B, by = k, suffixes = c(".l", ".b"))
cat("rows", nrow(m), "| warn or lost differ",
    sum(m$warn.l != m$warn.b | m$lost.l != m$lost.b), "\n")
# the battery's statement: lane farther from glm()/lmer() than the base
hi <- m$batt.l > m$batt.b & !is.na(m$batt.l) & !is.na(m$batt.b)
cat("rows where the lane is farther from the battery's reference:",
    sum(hi), "\n")
s <- function(x) if (all(is.na(x))) "NA" else
  sprintf("%.2e", max(x, na.rm = TRUE))
by_model <- do.call(rbind, lapply(split(m, m$model), function(x) {
  h <- x$batt.l > x$batt.b
  data.frame(model = x$model[1], n = nrow(x), batt_higher = sum(h),
             exact_lane_max = s(x$exact.l), exact_base_max = s(x$exact.b),
             lane_closer_exact = sum(x$exact.l < x$exact.b, na.rm = TRUE),
             lane_farther_exact = sum(x$exact.l > x$exact.b, na.rm = TRUE),
             exact_na = sum(is.na(x$exact.l) | is.na(x$exact.b)),
             batt_ref_err_max = s(abs(x$batt.l - x$exact.l)))
}))
print(by_model, row.names = FALSE)
cat("\nthe rows the battery called higher, against the exact reference:\n")
h <- m[hi, c(k, "batt.l", "batt.b", "exact.l", "exact.b")]
cat("of", nrow(h), ": lane closer to exact", sum(h$exact.l < h$exact.b,
                                                 na.rm = TRUE),
    "| farther", sum(h$exact.l > h$exact.b, na.rm = TRUE),
    "| equal", sum(h$exact.l == h$exact.b, na.rm = TRUE),
    "| NA", sum(is.na(h$exact.l) | is.na(h$exact.b)), "\n")
fmt <- h
for (v in c("batt.l", "batt.b", "exact.l", "exact.b")) {
  fmt[[v]] <- sprintf("%.3e", fmt[[v]])
}
print(fmt, row.names = FALSE)
