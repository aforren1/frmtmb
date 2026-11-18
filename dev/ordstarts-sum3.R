# Lane ordstarts fix round: the corrected counts (review
# dev/reviews/2026-11-17-ordstarts.md, items 1, 3, 4), from the logs of
# the main run (dev/ordstarts-out/main-*.tsv), the realness check
# (real.tsv, pilot-refits.tsv) and dev/ordstarts-swap.R (swap-*.tsv, which
# holds the default fit's warnings at summary()).
#
# silent: the default fit A gave no warning at the fit AND none at its
#   first summary().
# clean: the better point B has code 0, no degenerate warning, no lost
#   standard error ("Standard errors are not available"), every
#   standard error of the forced covariance at most 100, and holds under
#   the jittered refits.
#   Rscript dev/ordstarts-sum3.R > dev/ordstarts-out/sum3.txt
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE, quote = "",
                                    comment.char = "")
O <- "dev/ordstarts-out/"
X <- do.call(rbind, lapply(Sys.glob(paste0(O, "main-*.tsv")), rd))
key <- function(D) paste(D$family, D$re, D$n, D$seed)
A <- X[X$arm == "default", ]
B <- X[X$arm == "starts10", ]
B <- B[match(key(A), key(B)), ]
S <- do.call(rbind, lapply(Sys.glob(paste0(O, "swap-cum*.tsv")), rd))
SA <- S[S$what == "A", ]
iS <- match(key(A), key(SA))
cat("default fits with a summary() read:", sum(!is.na(iS)), "of", nrow(A),
    "; of those the swap run's A reproduces the main run's logLik (1e-6):",
    sum(abs(SA$ll[iS] - A$ll) < 1e-6, na.rm = TRUE), "\n")
d <- B$ll - A$ll
up <- d > 1e-3 & !is.na(d)
fitwarnA <- A$nwarn > 0
sumwarnA <- SA$s_nwarn[iS] > 0
silent <- !fitwarnA & !(sumwarnA %in% TRUE)
lostB <- grepl("Standard errors are not available", B$warns)
sdzeroB <- grepl("Boundary [(]singular[)] fit", paste(B$warns, B$msgs))
# the jitter check and the largest standard error at B
R <- rbind(rd(paste0(O, "real.tsv")), rd(paste0(O, "pilot-refits.tsv")))
holds <- rep(NA, nrow(A))
bigse <- rep(NA_real_, nrow(A))
for (j in which(up)) {
  r <- R[key(R) == key(A)[j], ]
  f0 <- r$ll[r$check == "fit" & r$variant == "B"][1]
  jj <- r[r$check == "real" & r$variant == "B", ]
  ok <- jj$code %in% 0
  holds[j] <- !any(ok) || max(jj$ll[ok]) - f0 <= 1e-3
  nt <- r$note[r$check == "se" & r$variant == "B"][1]
  bigse[j] <- suppressWarnings(as.numeric(sub(" at .*", "",
    sub("^[^ ]+ se ", "", sub(";.*", "", nt)))))
}
clean <- up & B$code %in% 0 & !B$degen_warn & !lostB & bigse <= 100 &
  holds %in% TRUE
des <- paste(A$family, A$re, paste0("n", A$n))
cat("\n== Up moves and the corrected silent shortfall ==\n")
cat("up", sum(up), "; A silent at the fit", sum(up & !fitwarnA),
    "; silent at the fit and at summary()", sum(up & silent), "\n")
cat("up, A silent at fit but summary() warns:", sum(up & !fitwarnA &
    sumwarnA %in% TRUE), ":", paste(key(A)[up & !fitwarnA &
    sumwarnA %in% TRUE], collapse = "; "), "\n")
cat("B lost a standard error (up moves):", sum(up & lostB), ":",
    paste(key(A)[up & lostB], collapse = "; "), "\n")
ss <- clean & silent
cat("SILENT SHORTFALL (silent A, clean B):", sum(ss), "of", nrow(A),
    "; with (1 | g):", sum(ss & A$re == "int"), "of",
    sum(A$re == "int"), "\n")
cat("gain:", signif(range(d[ss]), 4), "median", signif(median(d[ss]), 4),
    "\n")
print(table(des[ss]))
cat("of those, B has a group sd at zero (boundary maximum):",
    sum(ss & sdzeroB), "; interior:", sum(ss & !sdzeroB), "\n")
cat("rows:\n")
for (j in which(ss)) {
  cat(sprintf("  %-22s d %+.4f %s\n", key(A)[j], d[j],
              if (sdzeroB[j]) "sd at zero" else "interior"))
}
cat("\n== Hidden threshold SE on deferred defaults (item 3) ==\n")
h <- A$se_deferred %in% TRUE & A$tse1 > 100 | A$se_deferred %in% TRUE &
  A$tse2 > 100
h <- h %in% TRUE
tseA <- suppressWarnings(pmax(A$tse1, A$tse2, na.rm = TRUE))
cat("deferred defaults with a forced threshold SE above 100:", sum(h),
    "; of those with the degenerate warning anyway:", sum(h & A$degen_warn),
    "; hidden:", sum(h & !A$degen_warn), ":",
    paste(key(A)[h & !A$degen_warn], collapse = "; "), "\n")
cat("hidden ones whose summary() warns:", sum(h & !A$degen_warn &
    sumwarnA %in% TRUE), "; their threshold SE:",
    paste(signif(tseA[h & !A$degen_warn], 3), collapse = ", "), "\n")
same <- abs(d) <= 1e-6 & !is.na(d)
cat("degenerate warning gained by the starts arm at the same logLik",
    "(|d| <= 1e-6; the starts' evaluations built the covariance):",
    sum(same & B$degen_warn & !A$degen_warn), ":",
    paste(key(A)[same & B$degen_warn & !A$degen_warn], collapse = "; "),
    "\n")
cat("\n== Fit time per size, (1 | g) designs ==\n")
ran <- !is.na(B$ms_starts)
for (nn in c(300, 1000)) {
  s <- A$re == "int" & A$n == nn
  cat(sprintf("n %4d: default median %.2f s (%d fits); starts median %.2f s over the %d fits where starts ran\n",
              nn, median(A$elapsed[s]), sum(s), median(B$elapsed[s & ran]),
              sum(s & ran)))
}
