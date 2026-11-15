# Lane silent2, item 5207: the lane's fits against the base's on the
# controls. Usage: Rscript dev/silent2-ordcmp.R set [set ...]
sets <- commandArgs(TRUE)
rd <- function(s, arm) {
  x <- utils::read.delim(sprintf("dev/silent2-out/ord-%s-%s.tsv", s, arm))
  x$set <- s
  x[, c("case", "set", "ll", "code", "degen", "se_t1", "se_t2",
        grep("^se_t3$", names(x), value = TRUE))]
}
B <- do.call(rbind, lapply(sets, rd, "base"))
L <- do.call(rbind, lapply(sets, rd, "lane"))
M <- merge(B, L, by = c("case", "set"), suffixes = c(".b", ".l"))
num <- function(v) suppressWarnings(as.numeric(v))
M$se.b <- apply(M[, grep("^se_t[123]\\.b$", names(M))], 1, function(r) {
  r <- num(r)
  if (all(is.na(r))) NA else max(r, na.rm = TRUE)
})
cat("sets", paste(sets, collapse = " "), "; fits", nrow(M), "of base",
    nrow(B), "lane", nrow(L), "\n")
mv <- abs(M$ll.l - M$ll.b) > 1e-6
cat("logLik moved by > 1e-6:", sum(mv, na.rm = TRUE), "\n")
if (any(mv, na.rm = TRUE)) {
  print(M[mv %in% TRUE, c("case", "ll.b", "ll.l", "code.b", "code.l",
                          "degen.b", "degen.l")], row.names = FALSE)
}
cat("code changed:", sum(M$code.b != M$code.l, na.rm = TRUE), "\n")
cat("degenerate warning: base", sum(M$degen.b), "lane", sum(M$degen.l),
    "; gained", sum(!M$degen.b & M$degen.l), "; lost",
    sum(M$degen.b & !M$degen.l), "\n")
g <- !M$degen.b & M$degen.l
cat("gained, with the base's largest threshold SE (fixef) per component:\n")
print(M[g, c("case", "ll.b", "ll.l", "code.b", "se.b",
             grep("^se_t[123]\\.b$", names(M), value = TRUE))][
  order(M$se.b[g]), ], row.names = FALSE)
anyna <- apply(M[, grep("^se_t[123]\\.b$", names(M))], 1, function(r) {
  any(is.nan(num(r)))
})
cat("gained with every component's SE finite and at most 100:",
    sum(g & !anyna & (M$se.b <= 100) %in% TRUE), "\n")
cat("lost:\n")
print(M[M$degen.b & !M$degen.l, c("case", "ll.b", "ll.l", "se.b")],
      row.names = FALSE)
q <- !M$degen.l & M$code.l == 0
cat("passed by the lane (code 0, no warning):", sum(q),
    "; their largest threshold SE on the base arm:",
    format(max(M$se.b[q], na.rm = TRUE), digits = 4), "\n")
