# Which spelling brms-shaped consumers use to read a family's link
# fields: `$link` or `[["link"]]`. Counts occurrences in every function
# body of the installed packages that read a foreign model's family.
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pat_dollar <- "[$]link(fun|inv|_[a-z]+)?\\b"
pat_brack <- "\\[\\[\"link(fun|inv|_[a-z]+)?\"\\]\\]"
for (pkg in c("insight", "emmeans", "marginaleffects", "brms")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, "not installed\n")
    next
  }
  ns <- asNamespace(pkg)
  n_dollar <- 0L
  n_brack <- 0L
  for (nm in ls(ns, all.names = TRUE)) {
    f <- get(nm, ns)
    if (!is.function(f)) next
    txt <- paste(deparse(body(f)), collapse = "\n")
    n_dollar <- n_dollar +
      sum(lengths(regmatches(txt, gregexpr(pat_dollar, txt))))
    n_brack <- n_brack +
      sum(lengths(regmatches(txt, gregexpr(pat_brack, txt))))
  }
  cat(sprintf("%-16s $link*: %4d   [[\"link*\"]]: %4d\n", pkg, n_dollar,
              n_brack))
}
