# Lane predraw2, backlog 4932: the side law's cliff at draw_edge_c. On
# the short check's designs with a known predictive law (R: (1 | g), 6
# groups of 20, group sd 0.1 and 0.3; S: (1 + x | g), 15 groups of 8,
# slope sd 0.1 and 0.2; a new level), and on the review's cliff family
# (lambda 0.5 to 2), the exact coverage of predict()'s 90% and 98%
# intervals under four maps of a narrowed side, on the same deviates:
#   cur      0.82.0's map: min(z, z r, knots) where r = t2 / (2 csd) <= c
#   normal   no narrowing (the normal law, 0.81.0's draws)
#   blend0   (1 - a) z + a min(z, z r, knots), a = 1 - r / c
#   blend25  the same with a = 1 for r <= 0.25, falling to 0 at r = c
#   blend40  the same from r = 0.4
# r is the side's own ratio (v[2] / 2): the distance in conditional sds at
# which the log-likelihood falls by 2, over 2. Only fits with a narrowed
# side are run under every map; the rest draw the same in all.
# Usage: Rscript dev/predraw2-cliff-cov.R <lib> <part> <out.tsv>
#   part: R, S or cliff
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
part <- args[2]
out <- args[3]
cur_map <- ns$draw_edge_map
blend_map <- function(r0) {
  force(r0)
  function(x, s, v) {
    if (length(v) == 1L) return(x * v)
    r <- v[2] / 2
    a <- min(1, max(0, (ns$draw_edge_c - r) / (ns$draw_edge_c - r0)))
    (1 - a) * x + a * cur_map(x, s, v)
  }
}
maps <- list(cur = cur_map, normal = function(x, s, v) x,
             blend0 = blend_map(0), blend25 = blend_map(0.25),
             blend40 = blend_map(0.4))
set_map <- function(f) {
  unlockBinding("draw_edge_map", ns)
  assign("draw_edge_map", f, envir = ns)
  lockBinding("draw_edge_map", ns)
}
cover <- function(p, m, s, lev) {
  a <- (1 - lev) / 2
  qq <- stats::quantile(p, c(a, 1 - a), names = FALSE)
  c(stats::pnorm((qq[2] - m) / s) - stats::pnorm((qq[1] - m) / s),
    qq[2] - qq[1])
}
fq <- function(e) suppressMessages(suppressWarnings(e))
one <- function(tag, fit, nd, m, st) {
  # the law predict() builds and caches
  fq(predict(fit, newdata = nd, ndraws = 2, allow_new_levels = TRUE))
  law <- fit$cache[["draw_edge_law"]]
  rs <- if (!is.list(law)) NA_real_ else {
    j <- law$narrowed
    vs <- c(law$lo[j], law$hi[j])
    vs <- vs[vapply(vs, length, 1L) > 1L]
    if (length(vs)) min(vapply(vs, function(v) v[2] / 2, 0)) else NA_real_
  }
  arms <- if (is.na(rs)) "cur" else names(maps)
  rows <- list()
  for (a in arms) {
    set_map(maps[[a]])
    if (is.environment(fit$cache)) fit$cache[["draw_check"]] <- NULL
    set.seed(99)
    p <- fq(predict(fit, newdata = nd, ndraws = 1000, summary = FALSE,
                    allow_new_levels = TRUE))
    c90 <- cover(p, m, st, 0.9)
    c98 <- cover(p, m, st, 0.98)
    rows[[a]] <- data.frame(tag = tag, r = rs, arm = a, cov90 = c90[1],
                            w90 = c90[2], cov98 = c98[1], w98 = c98[2])
  }
  set_map(cur_map)
  do.call(rbind, rows)
}
res <- list()
if (part == "R") {
  for (sdg in c(0.1, 0.3)) for (seed in 1:100) {
    set.seed(1000 + seed)
    d <- data.frame(g = factor(rep(1:6, each = 20)), x = stats::rnorm(120))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, sdg)[d$g] + stats::rnorm(120)
    fit <- tryCatch(fq(frm(bf(y ~ x + (1 | g)), data = d)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    res[[length(res) + 1L]] <- one(paste0("R", sdg, ":", seed), fit,
                                   data.frame(x = 0, g = factor("new")), 1,
                                   sqrt(sdg^2 + 1))
  }
} else if (part == "S") {
  for (sds in c(0.1, 0.2)) for (seed in 1:100) {
    set.seed(2000 + seed)
    G <- 15
    m <- 8
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    x = stats::rnorm(G * m))
    ug <- stats::rnorm(G, 0, 0.8)
    sg <- stats::rnorm(G, 0, sds)
    d$y <- 1 + 0.5 * d$x + ug[d$g] + sg[d$g] * d$x + stats::rnorm(G * m)
    fit <- tryCatch(fq(frm(bf(y ~ x + (1 + x | g)), data = d)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    res[[length(res) + 1L]] <- one(paste0("S", sds, ":", seed), fit,
                                   data.frame(x = 2, g = factor("new")), 2,
                                   sqrt(0.8^2 + 4 * sds^2 + 1))
  }
} else {
  set.seed(5)
  dd <- data.frame(g = factor(rep(1:20, each = 10)), x = stats::rnorm(200))
  a <- stats::rnorm(20, 0, 0.7)[dd$g]
  b <- stats::rnorm(20, 0, 0.05)[dd$g] * dd$x
  e <- stats::rnorm(200)
  for (lam in seq(0.5, 2, by = 0.1)) {
    dd$y <- 1 + 0.5 * dd$x + a + lam * b + e
    fit <- fq(frm(bf(y ~ x + (1 + x | g)), data = dd))
    # the family's own truth is not a known law; the reviewer's reference
    # width is about 5.70 (dev/reviews/ppdraw-short-cliff.R)
    res[[length(res) + 1L]] <- one(sprintf("cliff:%.1f", lam), fit,
                                   data.frame(x = 2, g = factor("new")),
                                   1, 1)
  }
}
res <- do.call(rbind, res)
utils::write.table(res, out, sep = "\t", row.names = FALSE, quote = FALSE)
for (i in seq_len(nrow(res))) {
  cat(sprintf("%-12s r %6.3f %-8s cov90 %.3f cov98 %.3f w90 %.3f w98 %.3f\n",
              res$tag[i], res$r[i], res$arm[i], res$cov90[i], res$cov98[i],
              res$w90[i], res$w98[i]))
}
