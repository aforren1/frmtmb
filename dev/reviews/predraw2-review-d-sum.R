# Review of lane predraw2: summary of predraw2-review-d1/d2 logs.
x <- grep("seed", c(readLines("dev/reviews/predraw2-review-d1.txt"),
                    readLines("dev/reviews/predraw2-review-d2.txt")), value = TRUE)
x <- unique(x)
num <- function(p) as.numeric(vapply(regmatches(x, regexec(p, x)), `[`, "", 2L))
d <- data.frame(design = sub(" .*", "", x),
  fire = grepl("fire TRUE", x),
  excess = num("excess +(-?[0-9.]+)"),
  uniq = num("unique draws +([0-9]+)"),
  cl = num("coef cov90 lane ([0-9.]+)"), cb = num("coef cov90 lane [0-9.]+ base ([0-9.]+)"),
  pl = num("pred cov90 lane ([0-9.]+)"), pb = num("pred cov90 lane [0-9.]+ base ([0-9.]+)"),
  warn = grepl("frmtmb_draw_pareto", x))
for (g in unique(d$design)) {
  s <- d[d$design == g, ]; f <- s[s$fire, ]
  cat(sprintf("%-20s fits %2d fired %2d (excess %s) warn %2d | fired: unique %4.0f coef cov90 lane %.3f base %.3f pred cov90 lane %.4f base %.4f | max unfired excess %.3f\n",
    g, nrow(s), nrow(f), if (nrow(f)) paste(range(f$excess), collapse = "-") else "-",
    sum(f$warn), mean(f$uniq), mean(f$cl), mean(f$cb), mean(f$pl), mean(f$pb),
    max(c(-Inf, s$excess[!s$fire]))))
}
