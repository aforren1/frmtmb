# Review of lane ordstarts: an independent summary of the lane's logs
# (dev/ordstarts-out/main-*.tsv, real.tsv, pilot-refits.tsv, paths.tsv,
# moves.tsv). Written without reading dev/ordstarts-sum.R or sum2.R.
#   Rscript dev/reviews/ordstarts-review-sum.R
O <- "dev/ordstarts-out/"
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE,
                                    quote = "", comment.char = "")
X <- do.call(rbind, lapply(Sys.glob(paste0(O, "main-*.tsv")), rd))
X$key <- paste(X$family, X$re, X$n, X$seed)
A <- X[X$arm == "default", ]
B <- X[X$arm == "starts10", ][match(A$key, X$key[X$arm == "starts10"]), ]
stopifnot(identical(A$key, B$key), nrow(A) == 320)
cat("fits per arm", nrow(A), "\n")
errA <- A$err != "" & !is.na(A$err)
errB <- B$err != "" & !is.na(B$err)
cat("errors: A", sum(errA), "B", sum(errB), "same rows",
    identical(errA, errB), "\n")
d <- B$ll - A$ll
moved <- abs(d) > 1e-6 & !is.na(d)
up <- d > 1e-3 & !is.na(d)
dn <- d < -1e-3 & !is.na(d)
cat("moved", sum(moved), "up", sum(up), "down", sum(dn),
    "moved but |d| <= 1e-3", sum(moved & !up & !dn), "\n")
silent <- up & A$nwarn == 0
cat("silent up (A no warning)", sum(silent),
    "; also A no message", sum(up & A$nwarn == 0 & A$nmsg == 0), "\n")
ran <- !is.na(B$ms_starts)
cat("ran", sum(ran), "reached", sum(B$ms_reached, na.rm = TRUE), "of",
    10 * sum(ran), "\n")
cat("code1A (code != 0, not NA)", sum(A$code %in% 1), "code NA",
    sum(is.na(A$code)), "\n")
by <- paste(A$family, A$re, A$n)
tab <- data.frame(
  n = tapply(rep(1, 320), by, sum),
  err = tapply(errA, by, sum),
  code1A = tapply(A$code %in% 1, by, sum),
  proper = tapply(A$cls %in% "proper", by, sum),
  step = tapply(A$cls %in% "step", by, sum),
  flat = tapply(A$cls %in% "flat", by, sum),
  ran = tapply(ran, by, sum),
  reached = tapply(B$ms_reached, by, sum, na.rm = TRUE),
  moved = tapply(moved, by, sum), up = tapply(up, by, sum),
  down = tapply(dn, by, sum), silent = tapply(silent, by, sum),
  degA = tapply(A$degen_warn %in% TRUE, by, sum),
  degB = tapply(B$degen_warn %in% TRUE, by, sum),
  hintA = tapply(A$re_hint %in% TRUE, by, sum),
  stepB = tapply(B$step_msg %in% TRUE, by, sum))
print(tab)
print(colSums(tab))
cat("\nclass A x B\n"); print(table(A$cls, B$cls, useNA = "ifany"))
cat("code A x B\n"); print(table(A$code, B$code, useNA = "ifany"))
cat("up d:", signif(range(d[up]), 4), "median", signif(median(d[up]), 4),
    "\ndown d:", signif(range(d[dn]), 4), "\n")
cat("int moved", sum(moved & A$re == "int"), "slope moved",
    sum(moved & A$re == "slope"), "\n")
cat("proper A -> lower B:", sum(dn & A$cls == "proper"), "; code 0 -> 1:",
    sum(A$code %in% 0 & B$code %in% 1), "\n")
cat("down moves: A class", paste(table(A$cls[dn]), names(table(A$cls[dn]))),
    "; degA with hint", sum(dn & A$degen_warn & A$re_hint), "\n")
cat("silent per design\n"); print(tapply(silent, by, sum))
cat("silent up moves where A is code 0:", sum(silent & A$code %in% 0),
    "; A class:", paste(names(table(A$cls[silent])), table(A$cls[silent])),
    "\n")

cat("\n== cost ==\n")
el <- data.frame(A = tapply(A$elapsed, by, sum), B = tapply(B$elapsed, by,
                                                           sum))
el$ratio <- el$B / el$A
r <- B$elapsed / A$elapsed
el$med_ran <- tapply(r[ran], by[ran], median)
el$min_ran <- tapply(r[ran], by[ran], min)
el$max_ran <- tapply(r[ran], by[ran], max)
el$ev_ran <- tapply((B$evals / A$evals)[ran], by[ran], median)
print(round(el, 2))
cat("ALL A", sum(A$elapsed), "B", sum(B$elapsed), "ratio",
    sum(B$elapsed) / sum(A$elapsed), "per-fit ran median", median(r[ran]),
    "range", range(r[ran]), "evals median", median((B$evals / A$evals)[ran],
                                                  na.rm = TRUE), "\n")
cat("per-fit ratio over ALL non-error fits (not only ran): median",
    median(r[!errA]), "\n")
cat("int design fit times: A median", median(A$elapsed[A$re == "int"]),
    "B (ran) median", median(B$elapsed[A$re == "int" & ran]), "\n")
for (nn in c(300, 1000)) {
  s <- A$re == "int" & A$n == nn
  cat("  int n", nn, ": A median", median(A$elapsed[s]), "B ran median",
      median(B$elapsed[s & ran]), "\n")
}
# a starts fit is a default fit plus the starts, so the extra cost is B - A
cat("starts' added time per ran fit, median s:",
    median((B$elapsed - A$elapsed)[ran]), "\n")

cat("\n== deferred SE check ==\n")
ok <- !errA
cat("deferred A", sum(A$se_deferred[ok] %in% TRUE), "of", sum(ok),
    "; B", sum(B$se_deferred[!errB] %in% TRUE), "of", sum(!errB), "\n")
tA <- suppressWarnings(pmax(A$tse1, A$tse2, na.rm = TRUE))
hid <- A$se_deferred %in% TRUE & tA > 100 & !(A$degen_warn %in% TRUE)
cat("deferred A with forced thres SE > 100 and no degen warning:",
    sum(hid, na.rm = TRUE), ":", paste(A$key[hid %in% TRUE],
                                         collapse = "; "), "\n")
cat("deferred A with forced thres SE > 100 (any):",
    sum(A$se_deferred %in% TRUE & tA > 100, na.rm = TRUE), "\n")

cat("\n== error rows ==\n")
print(A[errA, c("key", "err")])

# the moved list the lane wrote
M <- rd(paste0(O, "moves.tsv"))
mk <- paste(M$family, M$re, M$n, M$seed)
cat("\nmoves.tsv rows", nrow(M), "match my moved set:",
    setequal(mk, A$key[moved]), "\n")
cat("max |moves.tsv d - my d|", max(abs(M$d - d[match(mk, A$key)])), "\n")

cat("\n== section 3, my reading ==\n")
R <- rbind(rd(paste0(O, "real.tsv")), rd(paste0(O, "pilot-refits.tsv")))
R$key <- paste(R$family, R$re, R$n, R$seed)
cat("keys in real+pilot", length(unique(R$key)), "\n")
U <- data.frame(key = A$key[up], d = d[up], warnA = A$nwarn[up],
                codeB = B$code[up], degB = B$degen_warn[up],
                tseB = suppressWarnings(pmax(B$tse1, B$tse2,
                                             na.rm = TRUE))[up],
                maxseB = B$maxse[up], mineigB = B$mineig[up],
                maxgradB = B$maxgrad[up], mineigA = A$mineig[up],
                clsA = A$cls[up], clsB = B$cls[up],
                stringsAsFactors = FALSE)
U$jupB <- U$jupA <- NA
U$jnB <- NA
U$topB <- ""
for (j in seq_len(nrow(U))) {
  r <- R[R$key == U$key[j], ]
  for (w in c("A", "B")) {
    f0 <- r$ll[r$check == "fit" & r$variant == w][1]
    jj <- r[r$check == "real" & r$variant == w, ]
    okk <- jj$code %in% 0
    U[[paste0("jup", w)]][j] <- if (any(okk)) max(jj$ll[okk]) - f0 else NA
    if (w == "B") U$jnB[j] <- sum(okk)
  }
  U$topB[j] <- r$note[r$check == "se" & r$variant == "B"][1]
  fb <- r$ll[r$check == "fit" & r$variant == "B"][1]
  if (abs(fb - B$ll[up][j]) > 1e-6) cat("B refit differs", U$key[j], "\n")
}
U$top1 <- suppressWarnings(as.numeric(sub(" at .*", "", sub(
  "^[^ ]+ se ", "", sub(";.*", "", U$topB)))))
U$holdsB <- is.na(U$jupB) | U$jupB <= 1e-3
clean <- U$codeB %in% 0 & !U$degB & !(U$tseB > 100) & U$holdsB
interior <- clean & !(U$top1 > 100)
cat("up", nrow(U), "silent", sum(U$warnA == 0), "clean", sum(clean),
    "silent+clean", sum(clean & U$warnA == 0), "silent+clean+all SE<=100",
    sum(interior & U$warnA == 0), "\n")
s24 <- interior & U$warnA == 0
cat("gain range of those", signif(range(U$d[s24]), 4), "\n")
cat("per design:\n"); print(table(sub(" [0-9]+$", "", U$key[s24])))
cat("B warns degenerate from silent A:", paste(U$key[U$degB & U$warnA == 0],
                                              collapse = "; "), "\n")
cat("jitter: B with no converged jittered refit (holds by default):",
    sum(U$jnB == 0), "\n")
cat("of the 24: min eigenvalue of B's outer Hessian (nll) min",
    signif(min(U$mineigB[s24], na.rm = TRUE), 3), "NA", sum(is.na(
      U$mineigB[s24])), "; max gradient max",
    signif(max(U$maxgradB[s24], na.rm = TRUE), 3), "\n")
cat("of the 24: A's min eigenvalue min", signif(min(U$mineigA[s24],
                                                   na.rm = TRUE), 3),
    "; A class", paste(names(table(U$clsA[s24])), table(U$clsA[s24])),
    "\n")
print(U[s24, c("key", "d", "tseB", "top1", "mineigB", "maxgradB",
               "mineigA", "jupA", "jupB", "jnB")], row.names = FALSE)
cat("\nUp moves not in the 24:\n")
print(U[!s24, c("key", "d", "warnA", "codeB", "degB", "tseB", "top1",
                "jupB", "jnB")], row.names = FALSE)
