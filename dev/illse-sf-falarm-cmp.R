# Lane illse: copy of dev/rel070-seflags-falarm-cmp.R. Row-by-row comparison of two false-alarm battery runs
# (dev/seflags-rev3-falarm.R's TSVs) on code, lost set, SE warnings,
# other warnings and max |SE/ref - 1|.
#   Rscript rel070-seflags-falarm-cmp.R <dirA> <tagA> <dirB> <tagB>
# reads <dir>/<group>-<cfg>-<tag>.tsv for cfg in ref, ob.
a <- commandArgs(trailingOnly = TRUE)
rd <- function(f) utils::read.delim(f, colClasses = "character")
for (cfg in strsplit(Sys.getenv("ILLSE_CFGS", "ref ob"), " ")[[1]]) {
  cat("\n==== BLAS", cfg, ":", a[2], "vs", a[4], "\n")
  tot <- 0; nd <- c(code = 0, lost = 0, se_warn = 0, other_warn = 0,
                    max_rel = 0)
  for (g in c("nl", "lin", "poly", "mixed", "bounded")) {
    x <- rd(file.path(a[1], sprintf("%s-%s-%s.tsv", g, cfg, a[2])))
    y <- rd(file.path(a[3], sprintf("%s-%s-%s.tsv", g, cfg, a[4])))
    stopifnot(identical(x$model, y$model), identical(x$seed, y$seed))
    tot <- tot + nrow(x)
    for (k in seq_len(nrow(x))) {
      d <- character()
      for (col in c("code", "lost", "se_warn", "other_warn")) {
        if (!identical(x[[col]][k], y[[col]][k])) {
          d <- c(d, col); nd[col] <- nd[col] + 1
        }
      }
      rx <- suppressWarnings(as.numeric(x$max_rel[k]))
      ry <- suppressWarnings(as.numeric(y$max_rel[k]))
      relchg <- !identical(x$max_rel[k], y$max_rel[k])
      if (relchg) nd["max_rel"] <- nd["max_rel"] + 1
      if (length(d) || relchg) {
        cat(sprintf(paste0("  %-20s seed %-3s | %-22s | code %s/%s | ",
                           "lost %s / %s | sew %s/%s | ow %s/%s | ",
                           "max_rel %s / %s\n"),
                    x$model[k], x$seed[k],
                    if (length(d)) paste(d, collapse = ",") else "max_rel only",
                    x$code[k], y$code[k], x$lost[k], y$lost[k],
                    x$se_warn[k], y$se_warn[k], x$other_warn[k],
                    y$other_warn[k], x$max_rel[k], y$max_rel[k]))
      }
    }
  }
  cat(sprintf("fits %d | rows differing in code %d, lost %d, se_warn %d, other_warn %d, max_rel (string) %d\n",
              tot, nd["code"], nd["lost"], nd["se_warn"], nd["other_warn"],
              nd["max_rel"]))
}
