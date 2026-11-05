# R3: the thresholds, by data set, from dev/latres-r3sum.R's table of fits.
x <- read.delim("dev/latres-out/r3-fits.tsv", stringsAsFactors = FALSE)
d <- do.call(rbind, lapply(split(x, x$design), function(z) data.frame(
  design = z$design[1], ratio = median(z$ratio, na.rm = TRUE), dar = median(z$dar, na.rm = TRUE),
  cov = mean(z$cov, na.rm = TRUE))))
d$bad <- d$ratio < .9 | d$ratio > 1.15 | abs(d$dar) > .15 | d$cov < .8
x$bad <- d$bad[match(x$design, d$design)]
few <- x$design %in% c("bern", "bern_long", "bern_big", "bern_s05", "ord4", "cum3", "cum5", "cum7", "binom3")
link_plain <- !x$design %in% c("pid_b", "pid_20", "pid_8", "bid_b", "ginv_b1", "ginv_b2", "ginv_b3", "ginv_s5",
                               "gident_b", "gid_s5", "ig_b1", "ig_b2", "ig_s5", "sn_u1", "sn_u10")
rule <- function(s0, sd, id) {
  ok <- x$g3 >= 1e-3
  w <- ok & (few | x$lambda < s0 | (x$disp > 0 & (x$lambda < sd | (link_plain & x$info_mean < id))))
  w[is.na(w)] <- FALSE
  w
}
res <- list()
for (s0 in c(.15, .18, .2)) for (sd in c(.25, .3, .35)) for (id in c(2, 2.5, 3)) {
  w <- rule(s0, sd, id)
  a <- x$set %in% c("lane", "m1b")
  res[[length(res) + 1]] <- data.frame(s0 = s0, sd = sd, id = id,
    hit = mean(w[a & x$bad]), fa = mean(w[a & !x$bad]),
    hit_all = mean(w[x$bad]), fa_all = mean(w[!x$bad]),
    nb_m20 = mean(w[x$design == "nb_m20"]), gamma = mean(w[x$design == "gamma"]),
    nb_s1 = mean(w[x$design == "nb_s1"]), nb = mean(w[x$design == "nb"]),
    zip_hi = mean(w[x$design == "zip_hi"]), pois = mean(w[x$design == "pois"]),
    pois_lo = mean(w[x$design == "pois_lo"]), pois_vlo = mean(w[x$design == "pois_vlo"]),
    ig_b2 = mean(w[x$design == "ig_b2"]), ginv_b2 = mean(w[x$design == "ginv_b2"]))
}
res <- do.call(rbind, res)
res <- res[order(res$fa - res$hit), ]
print(format(head(res[res$s0 %in% c(.15, .18, .2) & res$id %in% c(2.5, 3) & res$sd %in% c(.25, .3, .35), ], 12), digits = 2), row.names = FALSE)
