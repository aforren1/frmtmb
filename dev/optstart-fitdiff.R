# Lane optstart: which fits of the suites changed their optimum, their
# convergence code or their warnings between two runs of
# dev/optstart-scan.sh (dev/optstart-profile.R logs). Fits are matched
# within each test file by key and occurrence. An optimum changed when
# the objectives differ by more than 1e-8 relative
# (1e-8 * max(1, |objective|)); below that the two runs reached the
# same optimum to the optimizer's tolerance. Warnings are compared with
# their digits masked, so a number printed inside a message that moved
# in its last digits is not a new warning.
#   Rscript dev/optstart-fitdiff.R <base run> <lane run> [out prefix]
a <- commandArgs(TRUE)
rd <- function(run) {
  fs <- Sys.glob(file.path("dev/optstart-log", paste0("prof-", run), "*.txt"))
  L <- unlist(lapply(fs, function(f) {
    x <- readLines(f, warn = FALSE)
    paste0(basename(f), "\t", x)
  }))
  L <- L[grepl("\tFIT\t", L)]
  p <- strsplit(L, "\t", fixed = TRUE)
  g <- function(z, i) if (length(z) >= i) z[i] else NA_character_
  X <- data.frame(pid = vapply(p, g, "", 1), tag = vapply(p, g, "", 3),
                  key = vapply(p, g, "", 4), st = vapply(p, g, "", 5),
                  obj = suppressWarnings(as.numeric(vapply(p, g, "", 6))),
                  code = vapply(p, g, "", 7),
                  warn = vapply(p, g, "", 10),
                  form = vapply(p, g, "", 12),
                  stringsAsFactors = FALSE)
  X$warn[X$st != "OK"] <- NA
  X <- X[X$tag != "?", ]
  X$occ <- stats::ave(seq_len(nrow(X)), X$tag, X$key, FUN = seq_along)
  X
}
B <- rd(a[1])
N <- rd(a[2])
cat("fits logged:", a[1], nrow(B), "|", a[2], nrow(N), "\n")
cat("test files with a fit:", length(unique(B$tag)), "|",
    length(unique(N$tag)), "\n")
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
cat("\nsame optimum, convergence code changed:\n")
cc <- M$class == "same optimum" & M$codes %in% c("lane code != 0 only",
                                                 "base code != 0 only")
print(table(M$codes[cc]))
ch <- M[chg | cc | M$class %in% c("error base only", "error lane only"), ]
ch <- ch[order(ch$class, -abs(ch$d)), ]
cat("\nchanged fits (objective lane - base; negative is better):\n")
for (i in seq_len(nrow(ch))) {
  cat(sprintf("%-20s %-44s %+.4g (base %.10g code %s | lane %.10g code %s) %s\n",
              ch$class[i], ch$tag[i], ch$d[i], ch$obj.b[i], ch$code.b[i],
              ch$obj.n[i], ch$code.n[i], substr(ch$form.n[i] %||% "", 1, 70)))
}
# warnings, digits masked
mask <- function(w) {
  if (is.na(w) || w == "-") return(character(0))
  unique(gsub("[0-9]+([.][0-9]+)?(e[-+]?[0-9]+)?", "#",
              strsplit(w, " || ", fixed = TRUE)[[1]]))
}
cat("\nwarnings that appear or disappear (fits matched in both runs):\n")
nw <- 0L
for (i in which(ok2)) {
  wb <- mask(M$warn.b[i])
  wn <- mask(M$warn.n[i])
  gone <- setdiff(wb, wn)
  new <- setdiff(wn, wb)
  if (!length(gone) && !length(new)) next
  nw <- nw + 1L
  cat(sprintf("%-44s %s d=%+.3g codes %s/%s %s\n", M$tag[i], M$class[i], d[i],
              M$code.b[i], M$code.n[i], substr(M$form.n[i], 1, 60)))
  for (x in new) cat("    + ", substr(x, 1, 90), "\n", sep = "")
  for (x in gone) cat("    - ", substr(x, 1, 90), "\n", sep = "")
}
cat("fits whose warnings changed:", nw, "\n")
w <- M[M$class == "lane higher (worse)", ]
cat("worse by more than 1e-6 relative:",
    sum(w$d > 1e-6 * pmax(1, abs(w$obj.b))), "\n")
if (length(a) >= 3) {
  utils::write.table(M, paste0(a[3], ".tsv"), sep = "\t", quote = FALSE,
                     row.names = FALSE)
}
