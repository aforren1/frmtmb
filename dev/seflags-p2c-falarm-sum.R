# Lane seflags, punch round 2c: the reviewer's summary of the false-alarm
# battery (dev/seflags-rev3-falarm-sum.R) on dev/seflags-log/p2c/falarm,
# with the other warnings (the ill-conditioning one among them) counted too.
#   Rscript dev/seflags-p2c-falarm-sum.R
L <- "C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-log/p2c/falarm"
rd <- function(f) utils::read.delim(f, colClasses = "character")
for (cfg in c("ref", "ob")) {
  cat("\n==== BLAS", cfg, "\n")
  cat(sprintf("%-22s %4s | %-15s | %-15s | %-6s | %s\n", "model", "fits",
              "lost lane/base", "SEwarn l/b", "diffs",
              "max|SE/ref-1| lane (base) | ridge calls, min c0/ev"))
  for (g in c("nl", "lin", "poly", "mixed", "bounded")) {
    l <- rd(file.path(L, sprintf("%s-%s-lane.tsv", g, cfg)))
    b <- rd(file.path(L, sprintf("%s-%s-base.tsv", g, cfg)))
    stopifnot(identical(l$model, b$model), identical(l$seed, b$seed))
    for (m in unique(l$model)) {
      i <- l$model == m
      dl <- l$lost[i] != b$lost[i]
      ml <- suppressWarnings(max(as.numeric(l$max_rel[i]), na.rm = TRUE))
      mb <- suppressWarnings(max(as.numeric(b$max_rel[i]), na.rm = TRUE))
      nr <- sum(as.numeric(l$n_ridge[i]), na.rm = TRUE)
      mr <- suppressWarnings(min(as.numeric(l$min_c0_over_ev[i]),
                                 na.rm = TRUE))
      cat(sprintf("%-22s %4d | %6d / %-6d | %6d / %-6d | %-6d | %.3g (%.3g) | %d, %.3g\n",
                  m, sum(i), sum(l$lost[i] != "-"), sum(b$lost[i] != "-"),
                  sum(as.numeric(l$se_warn[i]), na.rm = TRUE),
                  sum(as.numeric(b$se_warn[i]), na.rm = TRUE), sum(dl),
                  ml, mb, nr, mr))
      ow <- which(l$other_warn[i] != b$other_warn[i])
      if (length(ow)) {
        cat(sprintf("    other warnings differ on %d fits (seeds %s)\n",
                    length(ow), paste(l$seed[i][ow], collapse = " ")))
      }
      for (j in which(dl)) {
        k <- which(i)[j]
        cat(sprintf("    seed %s code %s/%s: lane %s | base %s | rel %s/%s\n",
                    l$seed[k], l$code[k], b$code[k], l$lost[k], b$lost[k],
                    l$max_rel[k], b$max_rel[k]))
      }
    }
  }
}
