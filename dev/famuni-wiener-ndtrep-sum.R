# Summarize dev/famuni-out/wiener-ndtrep-scalar.txt and -hier.txt
# (written by dev/famuni-wiener-ndtrep.R). Counts come from the logs.
# Output: dev/famuni-out/wiener-ndtrep-sum.txt
dir <- "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-out"
rd <- function(f) {
  l <- readLines(file.path(dir, f))
  l <- l[grepl("^[0-9]{8} ", l)]
  l
}
cat("== scalar arm\n")
l <- rd("wiener-ndtrep-scalar.txt")
p <- strsplit(trimws(l), " +")
sc <- data.frame(seed = vapply(p, `[`, "", 1), n = vapply(p, `[`, "", 2),
                 arm = vapply(p, `[`, "", 3), conv = vapply(p, `[`, "", 4),
                 warn = as.integer(vapply(p, `[`, "", 5)),
                 ll = suppressWarnings(as.numeric(vapply(p, `[`, "", 6))),
                 ndt = suppressWarnings(as.numeric(vapply(p, `[`, "", 7))),
                 se = suppressWarnings(as.numeric(vapply(p, `[`, "", 8))),
                 el = suppressWarnings(as.numeric(vapply(p, `[`, "", 9))))
cat("rows", nrow(sc), " distinct seeds", length(unique(sc$seed)), "\n")
print(table(n = sc$n, arm = sc$arm, conv = sc$conv))
print(aggregate(warn ~ n + arm, sc, function(w) sum(w > 0)))
for (n in unique(sc$n)) {
  a <- sc[sc$n == n & sc$arm == "default", ]
  b <- sc[sc$n == n & sc$arm == "log", ]
  m <- merge(a, b, by = "seed", suffixes = c("_def", "_log"))
  ok <- m$conv_def == "0" & m$conv_log == "0"
  dll <- m$ll_log - m$ll_def
  dz <- abs(m$ndt_log - m$ndt_def) / m$se_def
  cat(sprintf(paste0("n = %s: paired %d, both conv 0 on %d; logLik log - ",
                     "default: min %.3e max %.3e; |dndt|/SE max %.3e; ",
                     "rows where log is lower by > 1e-6: %d, higher: %d\n"),
              n, nrow(m), sum(ok), min(dll[ok]), max(dll[ok]), max(dz[ok]),
              sum(dll[ok] < -1e-6), sum(dll[ok] > 1e-6)))
  cat(sprintf("   median elapsed s: default %.2f log %.2f\n",
              median(m$el_def), median(m$el_log)))
  bad <- m[!ok | abs(dll) > 1e-6, ]
  if (nrow(bad)) print(bad)
}
f <- file.path(dir, "wiener-ndtrep-hier.txt")
if (file.exists(f)) {
  cat("\n== hier arm\n")
  l <- rd("wiener-ndtrep-hier.txt")
  p <- strsplit(trimws(l), " +")
  h <- data.frame(seed = vapply(p, `[`, "", 1), arm = vapply(p, `[`, "", 2),
                  conv = vapply(p, `[`, "", 3),
                  warn = as.integer(vapply(p, `[`, "", 4)),
                  ll = suppressWarnings(as.numeric(vapply(p, `[`, "", 5))),
                  rmse = suppressWarnings(as.numeric(vapply(p, `[`, "", 6))),
                  sd = suppressWarnings(as.numeric(vapply(p, `[`, "", 7))),
                  el = suppressWarnings(as.numeric(vapply(p, `[`, "", 8))),
                  wtxt = vapply(p, function(z) paste(z[-(1:8)],
                                                     collapse = " "), ""))
  cat("rows", nrow(h), " distinct seeds", length(unique(h$seed)), "\n")
  print(table(arm = h$arm, conv = h$conv))
  print(aggregate(warn ~ arm, h, function(w) sum(w > 0)))
  print(aggregate(cbind(rmse, sd, el) ~ arm, h, median))
  print(aggregate(rmse ~ arm, h, function(r) round(range(r), 2)))
  print(table(h$arm, substr(h$wtxt, 1, 40)))
}
if (file.exists(f)) {
  w <- reshape(h[c("seed", "arm", "rmse", "sd")], idvar = "seed",
               timevar = "arm", direction = "wide")
  cat("\npaired: log rmse below ndt_group on", sum(w$rmse.log < w$rmse.ndt_group),
      "of", nrow(w), "; below global on", sum(w$rmse.log < w$rmse.global),
      "of", nrow(w), "\n")
  cat("sd(log ndt | s) under log link: median", median(w$sd.log), "range",
      paste(range(w$sd.log), collapse = " to "), "(truth 0.12)\n")
}
