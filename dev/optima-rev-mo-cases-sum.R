# Reviewer of lane optima: pairs dev/optima-rev-mo-cases.R's two logs.
#   Rscript dev/optima-rev-mo-cases-sum.R
rd <- function(f) {
  x <- grep("^CASE .* ll ", readLines(f), value = TRUE)
  fld <- function(k) sub(paste0(".* ", k, " ([^ ]+).*"), "\\1", x)
  data.frame(key = sub("^CASE ([^ ]+) seed ([0-9]+).*", "\\1 \\2", x),
             ll = as.numeric(fld("ll")), code = as.integer(fld("code")),
             ev = as.numeric(fld("evals")),
             nse = as.integer(fld("nonfinite_se")),
             warn = as.integer(fld("warn")), t = as.numeric(fld("t")),
             search = fld("search"), det = grepl("det TRUE TRUE", x))
}
b <- rd("dev/optima-rev-out/mo-cases-base.txt")
l <- rd("dev/optima-rev-out/mo-cases-lane.txt")
m <- merge(b, l, by = "key", suffixes = c(".b", ".l"))
m$d <- m$ll.l - m$ll.b
cat("paired fits", nrow(m), "\n")
cat("lane better by > 1e-6:", sum(m$d > 1e-6), "; worse by < -1e-6:",
    sum(m$d < -1e-6), "; min d", format(min(m$d), digits = 3), "\n")
print(m[abs(m$d) > 1e-6, c("key", "d", "code.b", "code.l", "search.l")])
cat("codes != 0: base", sum(m$code.b != 0), "lane", sum(m$code.l != 0), "\n")
cat("fits with a non-finite SE: base", sum(m$nse.b > 0), "lane",
    sum(m$nse.l > 0), "; with a warning: base", sum(m$warn.b > 0), "lane",
    sum(m$warn.l > 0), "\n")
cat("determinism (same process, twice): base", sum(m$det.b), "lane",
    sum(m$det.l), "of", nrow(m), "\n")
gain <- suppressWarnings(as.numeric(sub(".*/", "", m$search.l)))
cat("search ran on", sum(m$search.l != "-"), "lane fits; gained > 1e-6 on",
    sum(gain > 1e-6, na.rm = TRUE), "\n")
cat("evaluations: base", sum(m$ev.b), "lane", sum(m$ev.l), "ratio",
    format(sum(m$ev.l) / sum(m$ev.b), digits = 3), "\n")
g <- sub(" .*", "", m$key)
r <- tapply(seq_len(nrow(m)), g, function(i) sum(m$ev.l[i]) / sum(m$ev.b[i]))
print(round(sort(r), 2))
lg <- grepl("^large", m$key)
cat("large model: evals", paste(m$ev.b[lg], m$ev.l[lg], sep = "->"),
    "; seconds", paste(m$t.b[lg], m$t.l[lg], sep = "->"), "\n")
