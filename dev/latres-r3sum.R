# R3: the condition's candidates against what each design reads. Joins
# dev/latres-r3.R's statistics with each set's exact MLE. Writes r3-fits.tsv and
# prints one row per design.
rd <- function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"",
                             stringsAsFactors = FALSE, fill = TRUE, comment.char = "")
st <- do.call(rbind, lapply(list.files("dev/latres-out/r3", "tsv$", full.names = TRUE), rd))
ex <- list()
truth <- c(bern = 1, bern_long = 1, binom3 = .8, pois = .5, pois_ar9 = .3, pois_s1 = 1, pois_big = .5,
           bern_s05 = .5, bern_big = 1, binom5 = .8, binom10 = .8, binom20 = .8, binom10_lo = .8,
           pois_lo = .5, pois_vlo = .8, pois_3 = .5, nb = .5, nb_lo = .5, gamma = .5, ord4 = 1,
           zip = .5, zip_hi = .5, zip_lo = .5, hp = .5, zib = .6, cum3 = 1, cum5 = 1, cum7 = 1,
           nb_m1 = .5, nb_m5 = .5, nb_m20 = .5, nb_s1 = 1, beta_lo = .5, binom4 = .8, ln_s2 = 1,
           ln_s05 = .5, sn_u1 = .5, sn_u10 = 5, bb_T = .5, zoib_T = .5, glog_c = .1, glog_s5 = .15,
           glog_w = .05, blog_w = .2, ginv_b1 = .1, ginv_b2 = .2, ginv_b3 = .3, ginv_s5 = .3,
           gident_b = 2, gid_s5 = 2, ig_b1 = .3, ig_b2 = .5, ig_s5 = .5, pid_b = 4, pid_20 = 2,
           pid_8 = 1.2, bid_b = .04)
a <- list.files("dev/latres-out/m1", "tsv$", full.names = TRUE)
e1 <- do.call(rbind, lapply(a, function(f) {
  x <- rd(f)
  data.frame(design = x$design, seed = x$seed, x_sderr = pmin(x$x_sderr, 10), x_ar = x$x_ar)
}))
a <- list.files("dev/review-latres-out/m1b", "tsv$", full.names = TRUE)
e2 <- do.call(rbind, lapply(a, function(f) {
  x <- rd(f)
  data.frame(design = x$design, seed = x$seed, x_sderr = exp(x$x_ls), x_ar = x$x_ar)
}))
a <- list.files("dev/review-latres-out/m4b2", "tsv$", full.names = TRUE)
e3 <- do.call(rbind, lapply(a, function(f) {
  x <- rd(f)
  data.frame(design = x$design, seed = x$seed, x_sderr = exp(x$x_ls), x_ar = x$x_ar)
}))
x <- merge(st, rbind(e1, e2, e3), by = c("design", "seed"), all.x = TRUE)
x <- x[!is.na(x$sderr), ]
x$reason[is.na(x$reason) | x$reason == ""] <- "-"
x$ratio <- x$sderr / x$x_sderr
x$dar <- x$ar - x$x_ar
x$cov <- abs(log(x$sderr) - log(truth[x$design])) <= 1.96 * x$se_ls
x$read <- abs(log(x$ratio)) < log(1.25) & abs(x$dar) < .3
write.table(x, "dev/latres-out/r3-fits.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
tab <- do.call(rbind, lapply(split(x, x$design), function(z) data.frame(
  set = z$set[1], design = z$design[1], n = nrow(z), ratio = median(z$ratio, na.rm = TRUE),
  dar = median(z$dar, na.rm = TRUE), cov = mean(z$cov, na.rm = TRUE), read = mean(z$read, na.rm = TRUE),
  info = median(z$info_med), infom = median(z$info_mean), v = median(z$v), kap = median(z$kappa_med),
  lam = median(z$lambda), g3 = median(z$g3), warned = sum(z$warned, na.rm = TRUE),
  why = names(sort(table(z$reason), decreasing = TRUE))[1L], disp = z$disp[1])))
tab$bad <- tab$ratio < .9 | tab$ratio > 1.15 | abs(tab$dar) > .15 | tab$cov < .8
tab <- tab[order(tab$set, tab$bad, tab$ratio), ]
print(format(tab[, c("set","design","ratio","dar","cov","read","info","lam","g3","disp","bad","warned","why")], digits = 3), row.names = FALSE)
x$bad <- tab$bad[match(x$design, tab$design)]
for (st0 in list("lane", "m1b", c("lane", "m1b"), "m4b", c("lane", "m1b", "m4b"))) {
  z <- x[x$set %in% st0, ]
  cat(sprintf("%-16s hit %d of %d (%.3f); false alarm %d of %d (%.3f)\n", paste(st0, collapse = "+"),
              sum(z$warned & z$bad, na.rm = TRUE), sum(z$bad), mean(z$warned[z$bad], na.rm = TRUE),
              sum(z$warned & !z$bad, na.rm = TRUE), sum(!z$bad), mean(z$warned[!z$bad], na.rm = TRUE)))
}
