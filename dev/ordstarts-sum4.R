# Lane ordstarts fix round: the swapped start (dev/ordstarts-swap.R)
# against the 10 jittered starts of the main run, per design, from the
# logs swap-*.tsv, swapk3-*.tsv and main-*.tsv in dev/ordstarts-out.
#
# A good swap end: code 0, not a step function or at a flat boundary
# (mixture_ord_start_class()), no degenerate warning and no lost
# standard error at its fit.
#   Rscript dev/ordstarts-sum4.R > dev/ordstarts-out/sum4.txt
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE, quote = "",
                                    comment.char = "")
O <- "dev/ordstarts-out/"
key <- function(D) paste(D$family, D$re, D$n, D$seed)
X <- do.call(rbind, lapply(Sys.glob(paste0(O, "main-*.tsv")), rd))
A <- X[X$arm == "default", ]
B <- X[X$arm == "starts10", ]
B <- B[match(key(A), key(B)), ]
S <- do.call(rbind, lapply(Sys.glob(paste0(O, "swap-cum*.tsv")), rd))
SA <- S[S$what == "A", ]
good <- function(D) D$code %in% 0 & D$cls %in% "proper" & !(D$degW %in%
  TRUE) & !(D$lostW %in% TRUE)
fmt <- function(v) format(signif(v, 3))
cat("== K = 2: one swapped start from the default fit ==\n")
cat("up: the starts' B above A by > 1e-3. reachB: a good swap end at B's",
    "logLik or above (within 1e-3).\nswapUp: a good swap end above A by",
    "> 1e-3; aboveB: above B by > 1e-3. cost: swap refit elapsed over A's",
    "elapsed (same process), median.\n\n")
des <- unique(paste(A$family, A$re, A$n))
for (v in c("value", "raw", "either")) {
  cat("variant", v, "\n")
  for (dd in des) {
    a <- A[paste(A$family, A$re, A$n) == dd, ]
    b <- B[match(key(a), key(B)), ]
    sa <- SA[match(key(a), key(SA)), ]
    pick <- function(vv) {
      s <- S[S$what == "perm21" & S$variant == vv, ]
      s[match(key(a), key(s)), ]
    }
    if (v == "either") {
      s1 <- pick("value")
      s2 <- pick("raw")
      g1 <- good(s1)
      g2 <- good(s2)
      l1 <- ifelse(g1, s1$ll, -Inf)
      l2 <- ifelse(g2, s2$ll, -Inf)
      sll <- pmax(l1, l2)
      sgood <- is.finite(sll)
      el <- s1$elapsed + s2$elapsed
    } else {
      s <- pick(v)
      sgood <- good(s)
      sll <- ifelse(sgood, s$ll, -Inf)
      el <- s$elapsed
    }
    up <- (b$ll - a$ll > 1e-3) %in% TRUE
    reach <- sgood & sll >= b$ll - 1e-3
    sup <- sgood & sll - a$ll > 1e-3
    aboveB <- sgood & sll - b$ll > 1e-3
    cost <- el / sa$elapsed
    cat(sprintf(paste("  %-20s seeds %2d | up %2d reachB on up %2d | swapUp",
                      "%2d (where B not up %d) | aboveB %d | cost median %s",
                      "(range %s to %s)\n"),
                dd, sum(!is.na(a$ll)), sum(up), sum(reach & up, na.rm = TRUE),
                sum(sup, na.rm = TRUE), sum(sup & !up, na.rm = TRUE),
                sum(aboveB, na.rm = TRUE),
                fmt(median(cost, na.rm = TRUE)),
                fmt(min(cost, na.rm = TRUE)), fmt(max(cost, na.rm = TRUE))))
  }
}
cat("\nsame family (cum_cum): the permuted start's objective less A's (dfn)",
    "on code-0 defaults, and the swap end's logLik less A's:\n")
cc <- S[S$family == "cum_cum" & S$what == "perm21", ]
ccA <- SA[match(key(cc), key(SA)), ]
# the value variant gives each gap at least 0.1, so it is not the same
# point when a gap is smaller; the raw variant is the exact relabeling
for (v in c("raw", "value")) {
  o <- ccA$code %in% 0 & cc$variant == v
  cat("  ", v, ": refits", sum(o), "; |dfn| max",
      fmt(max(abs(cc$dfn[o]), na.rm = TRUE)), "median",
      fmt(median(abs(cc$dfn[o]), na.rm = TRUE)), "; dfn NA",
      sum(is.na(cc$dfn[o])), "; |ll - llA| max",
      fmt(max(abs(cc$ll[o] - ccA$ll[o]), na.rm = TRUE)), "\n")
}
cat("cum_sratio (different families): dfn on code-0 defaults, median",
    {
      cs <- S[S$family == "cum_sratio" & S$what == "perm21", ]
      csA <- SA[match(key(cs), key(SA)), ]
      o2 <- csA$code %in% 0
      paste(fmt(median(cs$dfn[o2], na.rm = TRUE)), "range",
            fmt(min(cs$dfn[o2], na.rm = TRUE)), "to",
            fmt(max(cs$dfn[o2], na.rm = TRUE)))
    }, "\n")

cat("\n== K = 3 (k3: cumulative, sratio, acat; (1 | g), n 1000, G 50;",
    "value variant, 5 permutations) ==\n")
K3 <- do.call(rbind, lapply(Sys.glob(paste0(O, "swapk3-*.tsv")), rd))
if (!is.null(K3) && nrow(K3)) {
  k3A <- K3[K3$what == "A", ]
  cat("seeds", nrow(k3A), "; errors", sum(k3A$err != "" &
      !is.na(k3A$err)), "; A code 0", sum(k3A$code %in% 0), "; A class",
      paste(names(table(k3A$cls)), table(k3A$cls), collapse = ", "), "\n")
  best <- gain <- cst <- rep(NA_real_, nrow(k3A))
  nup <- rep(0L, nrow(k3A))
  for (j in seq_len(nrow(k3A))) {
    p <- K3[K3$seed == k3A$seed[j] & grepl("^perm", K3$what), ]
    if (!nrow(p) || is.na(k3A$ll[j])) next
    g <- good(p)
    nup[j] <- sum(g & p$ll - k3A$ll[j] > 1e-3)
    if (any(g)) gain[j] <- max(p$ll[g]) - k3A$ll[j]
    cst[j] <- sum(p$elapsed) / k3A$elapsed[j]
    cat(sprintf(paste("  seed %2d A %.4f code %s %s | good perms %d,",
                      "above A %d, best good gain %s | cost %.2f\n"),
                k3A$seed[j], k3A$ll[j], k3A$code[j], k3A$cls[j], sum(g),
                nup[j], if (is.na(gain[j])) "NA" else fmt(gain[j]),
                cst[j]))
  }
  cat("seeds where a permutation reaches a good point above A by > 1e-3:",
      sum(nup > 0), "of", sum(!is.na(cst)), "; cost (5 refits over the",
      "fit) median", fmt(median(cst, na.rm = TRUE)), "range",
      fmt(min(cst, na.rm = TRUE)), "to", fmt(max(cst, na.rm = TRUE)), "\n")
}

cat("\n== The 19 silent shortfalls of dev/ordstarts-sum3.R: does the swap",
    "reach B? ==\n")
s3 <- readLines(paste0(O, "sum3.txt"))
rows <- s3[grepl("^  cum_[a-z]+ [a-z]+ [0-9]+ [0-9]+ +d [+-]", s3)]
ks <- sub(" +d .*", "", sub("^  ", "", rows))
for (v in c("value", "raw")) {
  s <- S[S$what == "perm21" & S$variant == v, ]
  s <- s[match(ks, key(s)), ]
  b <- B[match(ks, key(B)), ]
  r <- good(s) & s$ll >= b$ll - 1e-3
  cat(" ", v, ": reaches B or above on", sum(r), "of", length(ks), "\n")
}
s1 <- S[S$what == "perm21" & S$variant == "value", ]
s1 <- s1[match(ks, key(s1)), ]
s2 <- S[S$what == "perm21" & S$variant == "raw", ]
s2 <- s2[match(ks, key(s2)), ]
b <- B[match(ks, key(B)), ]
r <- (good(s1) & s1$ll >= b$ll - 1e-3) | (good(s2) & s2$ll >= b$ll - 1e-3)
cat("  either: reaches B or above on", sum(r), "of", length(ks),
    "; missed:", paste(ks[!r], collapse = "; "), "\n")

cat("\n== Without random effects (re none; A = starts 0, B = today's default",
    "of 10 starts; swapnore-*.tsv) ==\n")
NR <- do.call(rbind, lapply(Sys.glob(paste0(O, "swapnore-*.tsv")), rd))
if (!is.null(NR) && nrow(NR)) {
  for (dd in unique(paste(NR$family, NR$n))) {
    z <- NR[paste(NR$family, NR$n) == dd, ]
    a <- z[z$what == "A", ]
    b <- z[z$what == "B", ][match(a$seed, z$seed[z$what == "B"]), ]
    sv <- lapply(c(value = "value", raw = "raw"), function(v) {
      s <- z[z$what == "perm21" & z$variant == v, ]
      s[match(a$seed, s$seed), ]
    })
    up <- (b$ll - a$ll > 1e-3) %in% TRUE
    down <- (b$ll - a$ll < -1e-3) %in% TRUE
    gl <- lapply(sv, function(s) ifelse(good(s), s$ll, -Inf))
    sll <- pmax(gl$value, gl$raw)
    reach <- is.finite(sll) & sll >= b$ll - 1e-3
    sup <- is.finite(sll) & sll - a$ll > 1e-3
    aboveB <- is.finite(sll) & sll - b$ll > 1e-3
    cB <- b$elapsed / a$elapsed
    cR <- sv$raw$elapsed / a$elapsed
    cS <- (sv$raw$elapsed + sv$value$elapsed) / a$elapsed
    degNew <- (!a$degW & b$degW) %in% TRUE & up
    if (any(aboveB)) {
      cat("  ", dd, ": swap above B by", fmt(min(sll[aboveB] - b$ll[aboveB])),
          "to", fmt(max(sll[aboveB] - b$ll[aboveB])), "median",
          fmt(median(sll[aboveB] - b$ll[aboveB])), "
")
    }
    cat(sprintf(paste("  %-15s seeds %2d | B up %2d down %2d | B gains the",
                      "degenerate warning on an up move %d | either swap:",
                      "reachB on up %2d, swapUp %2d, aboveB %d | B/A sum",
                      "%.2f median %.2f | raw swap/A median %.2f | both",
                      "swaps/A median %.2f\n"),
                dd, nrow(a), sum(up), sum(down), sum(degNew),
                sum(reach & up), sum(sup), sum(aboveB),
                sum(b$elapsed) / sum(a$elapsed), median(cB, na.rm = TRUE),
                median(cR, na.rm = TRUE), median(cS, na.rm = TRUE)))
    for (v in c("value", "raw")) {
      g <- good(sv[[v]]) & sv[[v]]$ll >= b$ll - 1e-3
      cat(sprintf("      %-5s swap alone reaches B on up %d\n", v,
                  sum(g & up)))
    }
  }
}
