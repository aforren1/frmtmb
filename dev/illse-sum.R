# Lane illse: summarize dev/illse-check.R and dev/illse-glm.R logs, one
# row per design: fits, fits with a warning, fits with a lost SE, the
# largest entry error on the correlation scale and on the variances
# against the exact reference ("vs exact" for the glm logs), and the
# range of the fitted() SE ratio (to the exact SE, or to the poly() fit).
#   Rscript dev/illse-sum.R <log> [<log> ...]
num <- function(x) suppressWarnings(as.numeric(x))
for (f in commandArgs(TRUE)) {
  l <- readLines(f)
  l <- l[grepl(" s[0-9]+ code ", l)]
  if (!length(l)) next
  lab <- trimws(sub(" s[0-9]+ code .*", "", l))
  glm <- grepl("vs exact", l[1])
  cor <- num(sub(if (glm) ".*vs exact cor ([^ ]+).*" else
                   ".*cov cor ([^ ]+).*", "\\1", l))
  dg <- num(sub(if (glm) ".*vs exact cor [^ ]+ rel [^ ]+ diag ([^ ]+).*"
                else ".*diag ([^ ]+) \\|.*", "\\1", l))
  pc <- num(sub(".*vs poly cor ([^ ]+).*", "\\1", l))
  fr <- sub(".*fitted SE/(exact|poly) ", "", l)
  fr <- lapply(strsplit(trimws(fr), " +"), function(v) {
    num(v[grepl("^(NaN|NA|[-0-9.e+]+)$", v) & v != "-"])
  })
  warned <- !grepl("warn - ", l)
  lost <- !grepl("lost none", l)
  cat("\n==", f, "\n")
  cat(sprintf("%-32s %4s %5s %5s %9s %9s %9s %s\n", "design", "fits",
              "warn", "lost", "max cor", "max diag",
              if (glm) "poly cor" else "", "fitted SE ratio"))
  for (d in unique(lab)) {
    i <- which(lab == d)
    r <- unlist(fr[i])
    mx <- function(v) if (all(is.na(v))) NA else max(v, na.rm = TRUE)
    cat(sprintf("%-32s %4d %5d %5d %9.1e %9.1e %9s %s\n", d, length(i),
                sum(warned[i]), sum(lost[i]), mx(cor[i]), mx(dg[i]),
                if (glm) sprintf("%.1e", mx(pc[i])) else "",
                if (all(is.na(r))) "NaN" else
                  sprintf("%.5f to %.5f (%d NaN)", min(r, na.rm = TRUE),
                          max(r, na.rm = TRUE), sum(is.na(r)))))
  }
}
