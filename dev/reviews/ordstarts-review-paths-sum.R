# Review of lane ordstarts: section 4 (refit paths) recounted from
# dev/ordstarts-out/paths.tsv and pilot-refits.tsv.
#   Rscript dev/reviews/ordstarts-review-paths-sum.R
O <- "dev/ordstarts-out/"
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE,
                                    quote = "", comment.char = "")
P <- rd(paste0(O, "paths.tsv"))
L <- rd(paste0(O, "pilot-refits.tsv"))
key <- function(D) paste(D$family, D$re, D$n, D$seed)
P$key <- key(P)
L$key <- key(L)
cat("keys in paths:", length(unique(P$key)), "; pilot key:",
    unique(L$key), "\n")
cat("pilot key's checks in paths.tsv:",
    paste(unique(P$check[P$key == unique(L$key)]), collapse = ","), "\n")
# boot and refit from pilot for its key, influence from paths
Q <- rbind(L[L$check %in% c("fit", "boot", "refit", "allfit"), ],
           P[!(P$key %in% L$key & P$check %in% c("boot", "refit",
                                                 "allfit", "fit")), ])
pair <- function(chk, v1, v2) {
  x <- Q[Q$check == chk & Q$variant == v1 & Q$k > 0, ]
  y <- Q[Q$check == chk & Q$variant == v2 & Q$k > 0, ]
  m <- merge(x[, c("key", "k", "ll", "elapsed")],
             y[, c("key", "k", "ll", "elapsed")], by = c("key", "k"))
  dd <- m$ll.x - m$ll.y
  cat(sprintf("%-9s %-10s vs %-10s: pairs %d (both finite %d), above %d,",
              chk, v1, v2, nrow(m), sum(is.finite(dd)),
              sum(dd > 1e-3, na.rm = TRUE)),
      "below", sum(dd < -1e-3, na.rm = TRUE), "| elapsed",
      round(sum(x$elapsed), 1), round(sum(y$elapsed), 1), "\n")
  invisible(m)
}
pair("boot", "A_starts10", "A")
pair("boot", "B", "B_starts0")
pair("influence", "A_starts10", "A")
pair("influence", "B_starts0", "A")
pair("influence", "B_starts0", "A_starts10")
rf <- Q[Q$check == "refit", ]
cat("refit ms_starts by variant:\n")
print(tapply(rf$ms_starts, rf$variant, function(v) paste(unique(v),
                                                          collapse = ",")))
cat("refit elapsed sum by variant:\n")
print(round(tapply(rf$elapsed, rf$variant, sum), 1))
fits <- Q[Q$check == "fit", ]
for (kk in unique(rf$key)) {
  a <- fits$ll[fits$key == kk & fits$variant == "A"]
  b <- fits$ll[fits$key == kk & fits$variant == "B"]
  r <- rf[rf$key == kk, ]
  cat(kk, ": refit A - A", signif(r$ll[r$variant == "A"] - a, 3),
      "; refit B - B", signif(r$ll[r$variant == "B"] - b, 3),
      "; refit B_starts0 - B", signif(r$ll[r$variant == "B_starts0"] - b,
                                      3), "\n")
}
af <- Q[Q$check == "allfit", ]
cat("allfit rows", nrow(af[af$variant == "A", ]), "per A; any starts",
    any(af$ms_starts > 0, na.rm = TRUE), "\n")
for (kk in unique(af$key)) {
  b <- fits$ll[fits$key == kk & fits$variant == "B"]
  a <- fits$ll[fits$key == kk & fits$variant == "A"]
  x <- af[af$key == kk & af$variant == "A", ]
  cat(kk, ": allfit(A) best - A", signif(max(x$ll, na.rm = TRUE) - a, 3),
      "; reaches B?", any(abs(x$ll - b) < 1e-3, na.rm = TRUE), "\n")
}
# boot from B: within B's replicates, how often were starts run?
bb <- Q[Q$check == "boot" & Q$variant == "B", ]
cat("B boot replicates ms_starts:", paste(table(bb$ms_starts,
                                                useNA = "ifany"),
                                          names(table(bb$ms_starts,
                                                      useNA = "ifany")),
                                          collapse = "; "), "\n")
