# Lane ordstarts (backlog 6136): tables of the follow-up checks on the
# moved seeds, from the logs of dev/ordstarts-refits.R
# (dev/ordstarts-out/real.tsv, pilot-refits.tsv) and dev/ordstarts-paths.R
# (paths.tsv), joined with the main run (main-*.tsv).
#   Rscript dev/ordstarts-sum2.R > dev/ordstarts-out/sum2.txt
rd <- function(f) {
  if (!file.exists(f)) return(NULL)
  utils::read.delim(f, stringsAsFactors = FALSE, quote = "",
                    comment.char = "")
}
O <- "dev/ordstarts-out/"
X <- do.call(rbind, lapply(Sys.glob(paste0(O, "main-*.tsv")), rd))
R <- rbind(rd(paste0(O, "real.tsv")), rd(paste0(O, "pilot-refits.tsv")))
# the pilot's seed had its influence() run again with the third variant
# (dev/ordstarts-paths.R), which replaces the pilot's
pil <- rd(paste0(O, "pilot-refits.tsv"))
Pth <- rbind(pil[pil$check != "influence", ], rd(paste0(O, "paths.tsv")))
key <- function(D) paste(D$family, D$re, D$n, D$seed)
A <- X[X$arm == "default", ]
B <- X[X$arm == "starts10", ]
B <- B[match(key(A), key(B)), ]
tse <- function(D) suppressWarnings(pmax(D$tse1, D$tse2, na.rm = TRUE))
M <- utils::read.delim(paste0(O, "moves.tsv"), stringsAsFactors = FALSE)
k <- key(M)
iA <- match(k, key(A))
M$warnA <- A$nwarn[iA]
M$degB <- B$degen_warn[iA]
M$tseB <- tse(B)[iA]
M$llA <- A$ll[iA]
M$llB <- B$ll[iA]
fmt <- function(v, dg = 3) ifelse(is.na(v), "NA", format(signif(v, dg)))

cat("== Is each moved-to point a maximum? (dev/ordstarts-refits.R real) ==\n")
cat("jitter: 4 refits from the point, every outer coordinate moved by",
    "N(0, 0.1), no starts;\nup = largest code-0 end above the point",
    "(log-likelihood); se = the three largest standard errors\n\n")
M$holdsA <- M$holdsB <- NA
M$upA <- M$upB <- NA
M$seB <- M$seA <- ""
for (j in seq_len(nrow(M))) {
  r <- R[key(R) == k[j], ]
  if (!nrow(r)) next
  for (who in c("A", "B")) {
    f0 <- r$ll[r$check == "fit" & r$variant == who][1]
    jj <- r[r$check == "real" & r$variant == who, ]
    ok <- jj$code %in% 0
    up <- if (any(ok)) max(jj$ll[ok]) - f0 else NA
    M[[paste0("up", who)]][j] <- up
    M[[paste0("holds", who)]][j] <- is.na(up) || up <= 1e-3
    M[[paste0("se", who)]][j] <- r$note[r$check == "se" & r$variant == who][1]
  }
  # the refit reproduces the main run
  fb <- r$ll[r$check == "fit" & r$variant == "B"][1]
  if (abs(fb - M$llB[j]) > 1e-6) {
    cat("NOTE: B refit", fmt(fb, 10), "against main", fmt(M$llB[j], 10),
        "\n")
  }
}
M <- M[order(-(M$d > 0), M$family, M$re, M$n, M$seed), ]
for (j in seq_len(nrow(M))) {
  m <- M[j, ]
  cat(sprintf("%-10s %-5s n%-4d seed %2d d %+8.4f %s->%s | up A %s B %s\n",
              m$family, m$re, m$n, m$seed, m$d, m$clsA, m$clsB,
              fmt(m$upA), fmt(m$upB)))
  cat("    A:", m$seA, "\n    B:", m$seB, "\n")
}
cat("\nB holds under jitter:", sum(M$holdsB, na.rm = TRUE), "of",
    sum(!is.na(M$holdsB)), "; A holds:", sum(M$holdsA, na.rm = TRUE), "of",
    sum(!is.na(M$holdsA)), "\n")
# a variance at its zero boundary shows as a huge standard error of its
# log; that is a boundary maximum, not a flat direction of the mixture
first_name <- function(s) sub(" se .*", "", sub(";.*", "", s))
first_se <- function(s) suppressWarnings(as.numeric(
  sub(" at .*", "", sub("^[^ ]+ se ", "", sub(";.*", "", s)))))
M$bigB <- first_se(M$seB) > 100
M$bigB_name <- first_name(M$seB)
cat("\nB points whose largest standard error is above 100, by parameter:\n")
print(table(M$bigB_name[M$bigB %in% TRUE]))
U <- M[M$d > 1e-3, ]
cleanB <- U$codeB == 0 & !U$degB & !(U$tseB > 100) & U$holdsB %in% TRUE
sdB <- grepl("^theta|^sd_|log_sd|^theta_", U$bigB_name) & U$bigB %in% TRUE
cat("\n== Up moves (d > 1e-3) ==\n")
cat("up", nrow(U), "; silent default (no warning)", sum(U$warnA == 0),
    "\nB at code 0, no degenerate warning, threshold SE <= 100, holds",
    "under jitter:", sum(cleanB), "\n  of which the largest SE is above",
    "100:", sum(cleanB & U$bigB %in% TRUE), "(", paste(unique(
      U$bigB_name[cleanB & U$bigB %in% TRUE]), collapse = ", "), ")\n",
    "silent and clean as above:", sum(cleanB & U$warnA == 0),
    "; silent, clean, every SE <= 100:",
    sum(cleanB & U$warnA == 0 & !(U$bigB %in% TRUE)), "\n",
    "gain over those:", fmt(range(U$d[cleanB & U$warnA == 0 &
                                       !(U$bigB %in% TRUE)])), "\n")
U$kind <- ifelse(U$degB, "B warns degenerate",
  ifelse(!cleanB, "B other (code, SE or jitter)",
  ifelse(U$bigB %in% TRUE, "B clean but a log-gap or sd SE > 100",
         "B clean interior, every SE <= 100")))
U$kind[U$warnA != 0] <- "default warned"
cat("\nup moves by design and kind of the point B (defaults silent unless",
    "'default warned'):\n")
print(table(paste(U$family, U$re, paste0("n", U$n)), U$kind))
cat("\n")
cat("B with the degenerate warning (moved up from a silent default):",
    sum(U$degB & U$warnA == 0), ":", paste(with(U[U$degB & U$warnA == 0, ],
      sprintf("%s %s n%d seed %d thrSE %s", family, re, n, seed,
              fmt(tseB))), collapse = "; "), "\n")

cat("\n== Refit paths (dev/ordstarts-paths.R; cum_cum int n300 seed 27",
    "from dev/ordstarts-refits.R) ==\n")
pk <- unique(key(Pth))
for (kk in pk) {
  p <- Pth[key(Pth) == kk, ]
  llA <- p$ll[p$check == "fit" & p$variant == "A"][1]
  llB <- p$ll[p$check == "fit" & p$variant == "B"][1]
  cat("\n", kk, ": A", fmt(llA, 9), "B", fmt(llB, 9), "\n")
  rf <- p[p$check == "refit", ]
  for (v in unique(rf$variant)) {
    q <- rf[rf$variant == v, ]
    cat(sprintf("  refit %-10s ll %s starts run %s el %.1fs\n", v,
                fmt(q$ll, 9), q$ms_starts, q$elapsed))
  }
  bt <- p[p$check == "boot", ]
  if (nrow(bt)) {
    g <- function(v) bt$ll[bt$variant == v]
    e <- function(v) sum(bt$elapsed[bt$variant == v])
    dA <- g("A_starts10") - g("A")
    dB <- g("B") - g("B_starts0")
    cat(sprintf(paste("  boot A: replicates %d; A with starts above A by >",
                      "1e-3: %d (max %s), below: %d | el A %.1fs A10",
                      "%.1fs\n"),
                length(dA), sum(dA > 1e-3, na.rm = TRUE),
                fmt(suppressWarnings(max(dA, na.rm = TRUE))),
                sum(dA < -1e-3, na.rm = TRUE), e("A"), e("A_starts10")))
    cat(sprintf(paste("  boot B: B (starts inherited) above B with",
                      "starts 0 by > 1e-3: %d, below: %d, NA %d | starts run",
                      "in B's replicates %s | el B %.1fs B0 %.1fs\n"),
                sum(dB > 1e-3, na.rm = TRUE), sum(dB < -1e-3, na.rm = TRUE),
                sum(is.na(dB)), paste(unique(bt$ms_starts[bt$variant ==
                                                            "B"]),
                                      collapse = ","), e("B"),
                e("B_starts0")))
  }
  inf <- p[p$check == "influence" & p$k > 0, ]
  if (nrow(inf)) {
    g <- function(v) {
      x <- inf[inf$variant == v, ]
      x$ll[order(x$k)]
    }
    e <- function(v) sum(inf$elapsed[inf$variant == v])
    d1 <- g("A_starts10") - g("A")
    cat(sprintf(paste("  influence: deletions %d; A with starts above A by",
                      "> 1e-3: %d (max %s) | el A %.1fs A10 %.1fs\n"),
                length(d1), sum(d1 > 1e-3), fmt(max(d1)), e("A"),
                e("A_starts10")))
    if (length(g("B_starts0")) == length(d1)) {
      d2 <- g("B_starts0") - g("A_starts10")
      d3 <- g("B_starts0") - g("A")
      cat(sprintf(paste("  influence: B with starts 0 against A with starts:",
                        "above %d below %d (min %s); against A: above %d",
                        "| el B0 %.1fs\n"),
                  sum(d2 > 1e-3), sum(d2 < -1e-3), fmt(min(d2)),
                  sum(d3 > 1e-3), e("B_starts0")))
    }
  }
  af <- p[p$check == "allfit", ]
  if (nrow(af)) {
    for (v in c("A", "B")) {
      q <- af[af$variant == v, ]
      cat(sprintf(paste("  allfit %s: %d optimizers, ll range %s to %s,",
                        "any starts %s, any within 1e-3 of B %s\n"),
                  v, nrow(q), fmt(min(q$ll, na.rm = TRUE), 9),
                  fmt(max(q$ll, na.rm = TRUE), 9),
                  any(q$ms_starts > 0, na.rm = TRUE),
                  any(abs(q$ll - llB) < 1e-3, na.rm = TRUE)))
    }
  }
}

cat("\n== The degenerate warning's hint (starts = 10) on the main run ==\n")
hA <- A$re_hint %in% TRUE
cat("default fits with the hint:", sum(hA), "; with starts = 10:",
    "no degenerate warning", sum(hA & !B$degen_warn), "; starts ran",
    sum(hA & !is.na(B$ms_starts)), "; default code not 0",
    sum(hA & A$code != 0, na.rm = TRUE), "; default class step",
    sum(hA & A$cls == "step"), "flat", sum(hA & A$cls == "flat"), "\n")

cat("\n== Refit paths, totals over the seeds above ==\n")
tot <- function(check, v1, v0) {
  n <- above <- below <- 0
  for (kk in pk) {
    p <- Pth[key(Pth) == kk & Pth$check == check, ]
    if (check == "influence") p <- p[p$k > 0, ]
    x1 <- p[p$variant == v1, ]
    x0 <- p[p$variant == v0, ]
    if (!nrow(x1) || nrow(x1) != nrow(x0)) next
    d <- x1$ll[order(x1$k)] - x0$ll[order(x0$k)]
    n <- n + sum(!is.na(d))
    above <- above + sum(d > 1e-3, na.rm = TRUE)
    below <- below + sum(d < -1e-3, na.rm = TRUE)
  }
  el <- function(v) sum(Pth$elapsed[Pth$check == check & Pth$variant == v &
                                      (check != "influence" | Pth$k > 0)])
  cat(sprintf("%-9s %-10s vs %-10s: %d pairs, above by > 1e-3 %d, below %d | elapsed %s %.1fs, %s %.1fs\n",
              check, v1, v0, n, above, below, v1, el(v1), v0, el(v0)))
}
tot("boot", "A_starts10", "A")
tot("boot", "B", "B_starts0")
tot("influence", "A_starts10", "A")
tot("influence", "B_starts0", "A")
tot("influence", "B_starts0", "A_starts10")
rf <- Pth[Pth$check == "refit", ]
cat("refit: starts run on A", paste(unique(rf$ms_starts[rf$variant == "A"]),
    collapse = ","), "; on B", paste(unique(rf$ms_starts[rf$variant == "B"]),
    collapse = ","), "; B elapsed", sprintf("%.1f", sum(rf$elapsed[rf$variant ==
    "B"])), "s against B_starts0", sprintf("%.1f", sum(rf$elapsed[rf$variant ==
    "B_starts0"])), "s\n")
af <- Pth[Pth$check == "allfit", ]
cat("allfit: refits", nrow(af), "; any ran starts",
    any(af$ms_starts > 0, na.rm = TRUE), "\n")
