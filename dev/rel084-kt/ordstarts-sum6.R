# rel084-kt copy of dev/ordstarts-sum6.R: the after-*.tsv come from dev/rel084-log/kt/ordstarts-out/ (merged build); every before file is read from the ordstarts lane tree's dev/ordstarts-out/. No other change.
# Lane ordstarts, build round: before and after tables of C0 + D, from
# the logs (dev/ordstarts-out): the lane build's after-*.tsv
# (dev/ordstarts-after.R: A0 = starts 0, D = the new default), the main
# run main-*.tsv (rellib-r22: A = 0.83.0's default with random effects,
# B = starts 10), the realness check real.tsv / pilot-refits.tsv, the
# swap logs swap-*.tsv (summary() of A) and, without random effects,
# swapnore-*.tsv (rellib-r22: A = starts 0, B = 0.83.0's default of 10
# starts, the swapped ends).
#
# silent: no warning at the fit and none at the first summary().
# clean (with random effects, as dev/ordstarts-sum3.R): code 0, no
#   degenerate warning, no lost standard error, every standard error at
#   most 100, holds under the jitter check.
# clean (without random effects): code 0, proper, no degenerate and no
#   lost-standard-error warning (no jitter check there).
#   Rscript dev/ordstarts-sum6.R > dev/ordstarts-out/sum6.txt
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE, quote = "",
                                    comment.char = "")
O <- "C:/Users/adf44/source/r/frmtmb-wt-ordstarts/dev/ordstarts-out/"
key <- function(D) paste(D$family, D$re, D$n, D$seed)
fmt <- function(v) format(signif(v, 3))
AF <- do.call(rbind, lapply(Sys.glob("dev/rel084-log/kt/ordstarts-out/after-*.tsv"), rd))
A0 <- AF[AF$arm == "A0", ]
D <- AF[AF$arm == "D", ]
D <- D[match(key(A0), key(D)), ]
silentD <- D$nwarn == 0 & D$s_nwarn == 0
cat("== With random effects (8 designs, 40 seeds) ==\n")
X <- do.call(rbind, lapply(Sys.glob(paste0(O, "main-*.tsv")), rd))
A <- X[X$arm == "default", ]
B <- X[X$arm == "starts10", ]
B <- B[match(key(A), key(B)), ]
SA <- do.call(rbind, lapply(Sys.glob(paste0(O, "swap-cum*.tsv")), rd))
SA <- SA[SA$what == "A", ]
SA <- SA[match(key(A), key(SA)), ]
re <- A0$re != "none" & A0$family != "k3"
i <- match(key(A), key(A0))
cat("fits", nrow(A), "; lane A0 matches 0.83.0's default logLik (1e-9):",
    sum(abs(A0$ll[i] - A$ll) < 1e-9, na.rm = TRUE), "of",
    sum(!is.na(A$ll)), "; max |diff|", fmt(max(abs(A0$ll[i] - A$ll),
                                               na.rm = TRUE)), "\n")
d <- B$ll - A$ll
R <- rbind(rd(paste0(O, "real.tsv")), rd(paste0(O, "pilot-refits.tsv")))
cleanB <- rep(FALSE, nrow(A))
for (j in which(d > 1e-3)) {
  r <- R[key(R) == key(A)[j], ]
  f0 <- r$ll[r$check == "fit" & r$variant == "B"][1]
  jj <- r[r$check == "real" & r$variant == "B", ]
  ok <- jj$code %in% 0
  holds <- !any(ok) || max(jj$ll[ok]) - f0 <= 1e-3
  nt <- r$note[r$check == "se" & r$variant == "B"][1]
  big <- suppressWarnings(as.numeric(sub(" at .*", "",
    sub("^[^ ]+ se ", "", sub(";.*", "", nt)))))
  cleanB[j] <- B$code[j] %in% 0 & !B$degen_warn[j] &
    !grepl("Standard errors are not available", B$warns[j]) &
    big <= 100 & holds
}
silentA <- A$nwarn == 0 & !((SA$s_nwarn > 0) %in% TRUE)
Dm <- D[i, ]
A0m <- A0[i, ]
before <- silentA & cleanB & d > 1e-3
after <- (Dm$nwarn == 0 & Dm$s_nwarn == 0) & cleanB & (B$ll - Dm$ll > 1e-3)
des <- paste(A$family, A$re, paste0("n", A$n))
tab <- data.frame(design = unique(des))
for (k in seq_len(nrow(tab))) {
  s <- des == tab$design[k]
  tab$before[k] <- sum(before[s], na.rm = TRUE)
  tab$after[k] <- sum(after[s], na.rm = TRUE)
  tab$swap_kept[k] <- sum(Dm$swapped[s] == 1, na.rm = TRUE)
  tab$D_up[k] <- sum(Dm$ll[s] - A0m$ll[s] > 1e-3, na.rm = TRUE)
  tab$D_down[k] <- sum(Dm$ll[s] - A0m$ll[s] < -1e-6, na.rm = TRUE)
  tab$D_aboveB[k] <- sum(Dm$ll[s] - B$ll[s] > 1e-3, na.rm = TRUE)
  tab$degA0[k] <- sum(A0m$degW[s], na.rm = TRUE)
  tab$degD[k] <- sum(Dm$degW[s], na.rm = TRUE)
  r <- Dm$elapsed[s] / A0m$elapsed[s]
  tab$cost_sum[k] <- sum(Dm$elapsed[s], na.rm = TRUE) /
    sum(A0m$elapsed[s], na.rm = TRUE)
  tab$cost_med[k] <- stats::median(r, na.rm = TRUE)
  tab$cost_max[k] <- max(r, na.rm = TRUE)
}
print(tab, digits = 3, row.names = FALSE)
cat("ALL: silent shortfalls before", sum(before, na.rm = TRUE), "after",
    sum(after, na.rm = TRUE), "of", nrow(A), "; D down from a proper",
    "code-0 A0:", sum((A0m$cls == "proper" & A0m$code %in% 0 &
                         Dm$ll < A0m$ll - 1e-6) %in% TRUE),
    "; D down from any A0:", sum((Dm$ll < A0m$ll - 1e-6) %in% TRUE),
    "; cost D/A0 sum", fmt(sum(Dm$elapsed, na.rm = TRUE) /
                             sum(A0m$elapsed, na.rm = TRUE)),
    "median", fmt(stats::median(Dm$elapsed / A0m$elapsed, na.rm = TRUE)),
    "\n")
cat("after shortfalls:", paste(key(A)[after %in% TRUE], collapse = "; "),
    "\n")
cat("cost where the swap was kept: D/A0 median",
    fmt(stats::median((Dm$elapsed / A0m$elapsed)[Dm$swapped == 1],
                      na.rm = TRUE)),
    "; where a swap ran and was not kept:",
    fmt(stats::median((Dm$elapsed / A0m$elapsed)[Dm$swaps > 0 &
                                                   Dm$swapped == 0],
                      na.rm = TRUE)), "\n")

cat("\n== Without random effects (4 designs, 40 seeds) ==\n")
NR <- do.call(rbind, lapply(Sys.glob(paste0(O, "swapnore-*.tsv")), rd))
good <- function(Z) Z$code %in% 0 & Z$cls %in% "proper" &
  !(Z$degW %in% TRUE) & !(Z$lostW %in% TRUE)
for (dd in c("cum_cum 300", "cum_cum 1000", "cum_sratio 300",
             "cum_sratio 1000")) {
  z <- NR[paste(NR$family, NR$n) == dd, ]
  a <- z[z$what == "A", ]
  b <- z[z$what == "B", ][match(a$seed, z$seed[z$what == "B"]), ]
  sw <- z[grepl("^perm", z$what), ]
  dz <- D[D$re == "none" & paste(D$family, D$n) == dd, ]
  dz <- dz[match(a$seed, dz$seed), ]
  a0 <- A0[A0$re == "none" & paste(A0$family, A0$n) == dd, ]
  a0 <- a0[match(a$seed, a0$seed), ]
  best <- vapply(seq_len(nrow(a)), function(j) {
    cand <- c(if (good(b[j, ])) b$ll[j], if (good(dz[j, ])) dz$ll[j],
              sw$ll[sw$seed == a$seed[j] & good(sw)])
    if (length(cand)) max(cand) else NA_real_
  }, 0)
  silB <- b$nwarn == 0
  silD <- dz$nwarn == 0 & dz$s_nwarn == 0
  bef <- silB & (best - b$ll > 1e-3) %in% TRUE
  aft <- silD & (best - dz$ll > 1e-3) %in% TRUE
  cat(sprintf(paste("  %-16s | A0 = 0.83.0 starts 0 (1e-9): %d of %d |",
                    "short in silence: today's default %2d, D %2d | D vs",
                    "today: above %2d, below %2d (below and today clean %d)",
                    "| swap kept %2d | degenerate warning today %2d, D %2d",
                    "| cost today/A %.2f, D/A0 sum %.2f median %.2f\n"),
              dd, sum(abs(a0$ll - a$ll) < 1e-9, na.rm = TRUE), nrow(a),
              sum(bef), sum(aft), sum(dz$ll - b$ll > 1e-3, na.rm = TRUE),
              sum(dz$ll - b$ll < -1e-3, na.rm = TRUE),
              sum((dz$ll - b$ll < -1e-3) & good(b), na.rm = TRUE),
              sum(dz$swapped == 1, na.rm = TRUE), sum(b$degW, na.rm = TRUE),
              sum(dz$degW, na.rm = TRUE),
              sum(b$elapsed) / sum(a$elapsed),
              sum(dz$elapsed) / sum(a0$elapsed),
              stats::median(dz$elapsed / a0$elapsed)))
  cat("     D down from a proper code-0 A0:",
      sum((a0$cls == "proper" & a0$code %in% 0 & dz$ll < a0$ll - 1e-6) %in%
            TRUE), "\n")
}

cat("\n== K = 3 (cumulative, sratio, acat; n 3000; 20 seeds) ==\n")
k3 <- D[D$family == "k3", ]
k30 <- A0[A0$family == "k3", ]
k30 <- k30[match(k3$seed, k30$seed), ]
K3B <- rd(paste0(O, "swapk3s.tsv"))
k3b <- K3B[K3B$what == "B", ][match(k3$seed, K3B$seed[K3B$what == "B"]), ]
cat("swap kept", sum(k3$swapped == 1), "of", nrow(k3), "; D up over A0",
    sum(k3$ll - k30$ll > 1e-3), "; D down", sum(k3$ll - k30$ll < -1e-6),
    "; D above 10 starts (0.83.0)", sum(k3$ll - k3b$ll > 1e-3),
    "; below", sum(k3$ll - k3b$ll < -1e-3), "; cost D/A0 sum",
    fmt(sum(k3$elapsed) / sum(k30$elapsed)), "median",
    fmt(stats::median(k3$elapsed / k30$elapsed)), "\n")
cat("\nevery fit: errors A0", sum(A0$err != "" & !is.na(A0$err)), "D",
    sum(D$err != "" & !is.na(D$err)), "; step-start messages under D",
    sum(D$stepMsg %in% TRUE), "\n")
