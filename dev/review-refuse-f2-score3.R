# Joins R3 to the truth of review-refuse-f2-rule.R by tag and scores it beside
# the other rules. T0 must agree, which shows the data are the same.
O <- "dev/review-refuse-out/f2"
rd <- function(pat, key) {
  fs <- list.files(O, pattern = pat, full.names = TRUE)
  unlist(lapply(fs, function(f) grep(key, readLines(f, warn = FALSE), value = TRUE)))
}
a <- rd("^f2-rule-(table1|table2)", "^RULE t[12] .*truth_max")
b <- rd("^r3-(table1|table2)", "^R3 t[12] ")
g <- function(l, key) suppressWarnings(as.numeric(sub(paste0(".* ", key, " ([^ ]+).*"), "\\1", l)))
ta <- sub("^RULE (t[12] .*) truth_max.*", "\\1", a)
tb <- sub("^R3 (t[12] .*) T2 .*", "\\1", b)
d <- data.frame(tag = ta, truth = g(a, "truth_max"), T0a = g(a, "T0"), R1 = g(a, "R1"), stringsAsFactors = FALSE)
e <- data.frame(tag = tb, T2 = g(b, "T2"), mv = g(b, "mv"), warn = grepl("warn TRUE", b), T0b = g(b, "T0"), stringsAsFactors = FALSE)
m <- merge(d, e, by = "tag")
cat("joined", nrow(m), "of", nrow(d), "; T0 max |diff|", max(abs(m$T0a - m$T0b), na.rm = TRUE), "\n")
m$table <- substr(m$tag, 1, 2)
hit <- function(w, k) sprintf("hit %d of %d; false alarm %d of %d; warns on %d of %d at 0.5 to 1; largest shift silent %.2f; warns on %.0f%% of fits",
                              sum(k & m$truth >= 1 & w), sum(k & m$truth >= 1), sum(k & m$truth < 0.5 & w), sum(k & m$truth < 0.5),
                              sum(k & m$truth >= 0.5 & m$truth < 1 & w), sum(k & m$truth >= 0.5 & m$truth < 1), max(c(0, m$truth[k & !w])), 100 * mean(w[k]))
rules <- list(
  "current |T| > 1" = !is.na(m$T0a) & abs(m$T0a) > 1,
  "R1 gate+all own, 30 rows >= 0.3" = !is.na(m$T0a) & abs(m$T0a) > 1 & !is.na(m$R1) & m$R1 >= 0.3,
  "R3 TK gate, 10 rows, 0.3" = m$warn,
  "R3 at 0.5" = !is.na(m$mv) & m$mv >= 0.5,
  "R3 at 0.2" = !is.na(m$mv) & m$mv >= 0.2,
  "TK gate alone" = is.finite(m$T2) & abs(m$T2) > 1 | !is.finite(m$T2))
for (tb in c("t1", "t2", "all")) {
  k <- if (tb == "all") rep(TRUE, nrow(m)) else m$table == tb
  cat("---", tb, "\n")
  for (n in names(rules)) cat(sprintf("%-32s %s\n", n, hit(rules[[n]], k)))
}
cat("R3 misses at truth >= 0.75:\n"); print(m[m$truth >= 0.75 & !m$warn, c("tag", "truth", "T0a", "T2", "mv")])
cat("R3 false alarms below 0.3:\n"); print(m[m$truth < 0.3 & m$warn, c("tag", "truth", "T0a", "T2", "mv")])
cat(sprintf("R3 move vs truth, fits past the gate: cor %.3f; median |mv - truth| %.3f\n", cor(m$mv, m$truth, use = "complete.obs"), median(abs(m$mv - m$truth), na.rm = TRUE)))
# von_mises
v <- rd("^f2-rule-vm-x", "^VM mean")
r <- rd("^r3-vm-x", "^R3 vm")
if (length(v)) {
  tv <- sub("^VM mean ([0-9.]+) kappa ([0-9]+) read ([0-9]+) slope ([0-9.]+) rep ([0-9]).*", "\\1 \\2 \\3 \\4 \\5", v)
  vt <- data.frame(tag = tv, truth = suppressWarnings(as.numeric(sub(".*truth_max ([^;]+);.*", "\\1", v))), gap = as.numeric(sub(".*logLik - exact ([^;]+);.*", "\\1", v)),
                   warned = grepl("warned TRUE", v), R1 = suppressWarnings(as.numeric(sub(".*R1 ([^;]+);.*", "\\1", v))), pv = as.numeric(sub(".*my integrate ([^,]+),.*", "\\1", v)), stringsAsFactors = FALSE)
  tr <- sub("^R3 vm mean ([0-9.]+) kappa ([0-9]+) read ([0-9]+) slope ([0-9.]+) rep ([0-9]) T2.*", "\\1 \\2 \\3 \\4 \\5", r)
  rt <- data.frame(tag = tr, r3 = grepl("warn TRUE", r), mv = g(r, "mv"), T2 = g(r, "T2"), stringsAsFactors = FALSE)
  w <- merge(vt, rt, by = "tag", all.x = TRUE)
  w$end <- sub(" .*", "", w$tag) != "0.00"
  for (e in c(FALSE, TRUE)) {
    k <- w$end == e
    cat(sprintf("von_mises %s: %d fits; per value %.3f to %.3f; logLik - exact %.1f to %.1f; truth median %.2f max %.2f; truth >= 1 on %d; current warns %d; R3 warns %d of %d run\n",
                if (e) "near the end" else "mean 0", sum(k), min(w$pv[k]), max(w$pv[k]), min(w$gap[k]), max(w$gap[k]), median(w$truth[k], na.rm = TRUE), max(w$truth[k], na.rm = TRUE), sum(w$truth[k] >= 1, na.rm = TRUE), sum(w$warned[k]), sum(w$r3[k], na.rm = TRUE), sum(!is.na(w$r3[k]))))
  }
  print(w[w$truth >= 1 & !is.na(w$truth) & !isTRUE(w$r3), c("tag", "truth", "gap", "pv", "T2", "mv", "r3")][1:min(8, sum(w$truth >= 1, na.rm = TRUE)), ])
}
