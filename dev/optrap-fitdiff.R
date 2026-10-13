# Lane optrap: which fits of the suites changed their optimum between
# two runs of dev/optrap-scan.sh (dev/optrap-profile.R logs). Fits are
# matched within each test file by key and occurrence. A change is an
# objective difference above 1e-8 relative (1e-8 * max(1, |objective|));
# below it the two runs reached the same optimum to the optimizer's
# tolerance.
#   Rscript dev/optrap-fitdiff.R <base run> <lane run> [out prefix]
a <- commandArgs(TRUE)
rd <- function(run) {
  fs <- Sys.glob(file.path("dev/optrap-log", paste0("prof-", run), "*.txt"))
  L <- unlist(lapply(fs, function(f) {
    x <- readLines(f, warn = FALSE)
    paste0(basename(f), "\t", x)
  }))
  L <- L[grepl("\tFIT\t", L)]
  p <- strsplit(L, "\t", fixed = TRUE)
  X <- data.frame(pid = vapply(p, `[`, "", 1), tag = vapply(p, `[`, "", 3),
                  key = vapply(p, `[`, "", 4), st = vapply(p, `[`, "", 5),
                  obj = suppressWarnings(as.numeric(vapply(p, function(z)
                    if (length(z) >= 6) z[6] else NA_character_, ""))),
                  code = vapply(p, function(z)
                    if (length(z) >= 7) z[7] else NA_character_, ""),
                  msg = vapply(p, function(z) z[length(z)], ""),
                  stringsAsFactors = FALSE)
  X$line <- seq_len(nrow(X))
  X <- X[X$tag != "?", ]
  X$occ <- stats::ave(seq_len(nrow(X)), X$tag, X$key, FUN = seq_along)
  X
}
B <- rd(a[1])
N <- rd(a[2])
cat("fits logged:", a[1], nrow(B), "|", a[2], nrow(N), "\n")
cat("test files:", length(unique(B$tag)), "|", length(unique(N$tag)), "\n")
M <- merge(B, N, by = c("tag", "key", "occ"), all = TRUE,
           suffixes = c(".b", ".n"))
both <- !is.na(M$st.b) & !is.na(M$st.n)
ok2 <- both & M$st.b == "OK" & M$st.n == "OK" & is.finite(M$obj.b) &
  is.finite(M$obj.n)
d <- M$obj.n - M$obj.b
tol <- 1e-8 * pmax(1, abs(M$obj.b))
M$class <- ifelse(!both, ifelse(is.na(M$st.b), "lane only", "base only"),
  ifelse(M$st.b == "ERROR" & M$st.n == "ERROR", "error both",
    ifelse(M$st.b == "ERROR", "error base only",
      ifelse(M$st.n == "ERROR", "error lane only",
        ifelse(!ok2, "unfitted or non-finite",
          ifelse(abs(d) <= tol, "same optimum",
            ifelse(d < 0, "lane lower (better)", "lane higher (worse)")))))))
M$d <- d
# by convergence code: a change between two stops neither of which
# converged is a different stop, not a different maximum
cb <- suppressWarnings(as.integer(M$code.b))
cn <- suppressWarnings(as.integer(M$code.n))
M$codes <- ifelse(!ok2, "",
  ifelse(cb == 0 & cn == 0, "both code 0",
    ifelse(cb == 0, "lane code != 0 only",
      ifelse(cn == 0, "base code != 0 only", "neither code 0"))))
print(table(M$class))
cat("\nchanged optima by convergence codes:\n")
chg <- M$class %in% c("lane lower (better)", "lane higher (worse)")
print(table(M$class[chg], M$codes[chg]))
ch <- M[M$class %in% c("lane lower (better)", "lane higher (worse)",
                       "error base only", "error lane only"), ]
ch <- ch[order(ch$class, -abs(ch$d)), ]
cat("\nchanged fits (objective lane - base; negative is better):\n")
for (i in seq_len(nrow(ch))) {
  cat(sprintf("%-20s %-42s %+.4g (base %.10g code %s | lane %.10g code %s) %s\n",
              ch$class[i], ch$tag[i], ch$d[i], ch$obj.b[i], ch$code.b[i],
              ch$obj.n[i], ch$code.n[i],
              substr(if (is.null(ch$msg.n)) "" else ch$msg.n[i], 1, 70)))
}
w <- M[M$class == "lane higher (worse)", ]
cat("\nworse by more than 1e-6 relative:",
    sum(w$d > 1e-6 * pmax(1, abs(w$obj.b))), "\n")
if (length(a) >= 3) {
  utils::write.table(M, paste0(a[3], ".tsv"), sep = "\t", quote = FALSE,
                     row.names = FALSE)
}
