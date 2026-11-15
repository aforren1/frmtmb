# Summary of dev/predraw2-mixcov.R: lane against base per design, over the
# fits where the lane's check fired and over all fits.
# Usage: Rscript dev/predraw2-mixcov-sum.R <prefix> (reads
#        <prefix>-<g|p|c>-<lane|base>.tsv)
pre <- commandArgs(TRUE)[1]
for (d in c("g", "p", "c")) {
  fl <- paste0(pre, "-", d, "-lane.tsv")
  fb <- paste0(pre, "-", d, "-base.tsv")
  if (!file.exists(fl) || !file.exists(fb)) next
  a <- utils::read.delim(fl)
  b <- utils::read.delim(fb)
  m <- merge(a, b, by = "seed", suffixes = c(".l", ".b"))
  for (sel in c("fired", "not fired", "all")) {
    x <- switch(sel, fired = m[m$fire.l, ], "not fired" = m[!m$fire.l, ],
                all = m)
    if (!nrow(x)) next
    f <- function(v) sprintf("%.3f", mean(v))
    cat(sprintf(paste("%s %-9s n %2d | cov90 lane %s base %s plug %s |",
                      "cov98 lane %s base %s plug %s | mae lane %s base %s",
                      "plug %s | mean|z| lane %s base %s | k>0.7 %d\n"),
                d, sel, nrow(x), f(x$cov90.l), f(x$cov90.b),
                f(x$cov90_plug.l), f(x$cov98.l), f(x$cov98.b),
                f(x$cov98_plug.l), sprintf("%.4f", mean(x$mae_pred.l)),
                sprintf("%.4f", mean(x$mae_pred.b)),
                sprintf("%.4f", mean(x$mae_plug.l)), f(abs(x$z.l)),
                f(abs(x$z.b)),
                sum(suppressWarnings(as.numeric(sub("[.]$", "", x$k.l))) >
                      0.7 |
                      grepl("few", x$k.l), na.rm = TRUE)))
  }
  same <- m[!m$fire.l, ]
  cat(sprintf("%s unfired fits with summaries equal to 1e-9 lane/base: %d of %d\n",
              d, sum(abs(same$mae_pred.l - same$mae_pred.b) < 1e-9 &
                       (abs(same$z.l - same$z.b) < 1e-6 | is.na(same$z.l)), na.rm = TRUE),
              nrow(same)))
  cat(sprintf("%s excess of fired fits: %s\n", d,
              paste(sprintf("%.2f", sort(m$excess.l[m$fire.l])),
                    collapse = " ")))
}
