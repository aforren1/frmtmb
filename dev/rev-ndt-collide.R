# rev-ndt round 2, attack 1: does ddm_label_code() collide?
#
# The floor table is now keyed on a code derived from the group's LABEL:
# two polynomial rolling hashes, base 131 mod 33554393 and base 137 mod
# 16777213 started at 1, packed as h1 * 2^24 + h2. Distinct labels
# sharing a code inside one column are refused by ddm_coerce_ndt_group();
# a fitted label and a newdata-only label that collide are not, and that
# is the case this looks for.
#
# Three things are measured, in order of what a reader should believe:
#
#   1. a proof, not a search, that NO two labels of three or fewer
#      printable ASCII characters can collide;
#   2. a brute-force hunt over 10^8 labels of the shape a participant
#      identifier actually has, "S" plus eight digits, where the
#      birthday count says to expect about nine colliding pairs;
#   3. the collision found, checked against the package's own
#      ddm_label_code() rather than against this script's copy of it.
#
# No seed: the search is exhaustive over a stated set.

.libPaths(c("C:/Users/adf44/source/r/rev-ndt-lib2",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.eam))

P1 <- 33554393
P2 <- 16777213

# the package's own, reached through its namespace
pkg_code <- get("ddm_label_code", envir = asNamespace("frmtmb.eam"))

# the same function, vectorized over a matrix of code points, one row
# per label. Checked against the package's before it is trusted.
vec_code <- function(cpm) {
  h1 <- numeric(nrow(cpm))
  h2 <- rep(1, nrow(cpm))
  for (j in seq_len(ncol(cpm))) {
    k <- cpm[, j]
    h1 <- (h1 * 131 + k) %% P1
    h2 <- (h2 * 137 + k) %% P2
  }
  h1 * 16777216 + h2
}

chk <- c("s1", "subject-07", "Muller", "", "ZZZZ", "S00000000")
cpm <- t(vapply(chk, function(s) {
  cp <- utf8ToInt(s)
  if (!length(cp)) cp <- 0L
  c(cp, rep(NA_integer_, 10L - length(cp)))
}, numeric(10)))
mine <- vapply(seq_along(chk), function(i) {
  v <- cpm[i, ]
  vec_code(matrix(v[!is.na(v)], nrow = 1L))
}, numeric(1))
theirs <- pkg_code(chk)
cat("vectorized hash agrees with the package's:",
    identical(mine, theirs), "\n")
if (!identical(mine, theirs)) {
  print(data.frame(label = chk, mine = mine, theirs = theirs))
  stop("the reimplementation is wrong; nothing below is evidence")
}

# --------------------------------------------------------------------
# 1. Short ASCII labels cannot collide, and it is arithmetic, not luck.
#
# For a label of n printable ASCII characters (code points 32 to 126)
# h1 = sum cp_k 131^(n-k) and h2 = 137^n + sum cp_k 137^(n-k). If
# neither sum reaches its modulus there is no reduction at all, and a
# base-131 (or base-137) numeral with digits below the base is unique.
cat("\n-- 1. the no-wrap bound for printable ASCII\n")
for (n in 1:5) {
  m1 <- 126 * sum(131^(seq_len(n) - 1))
  m2 <- 137^n + 126 * sum(137^(seq_len(n) - 1))
  cat(sprintf("  n=%d  max h1 %.0f vs P1 %d %s |  max h2 %.0f vs P2 %d %s\n",
              n, m1, P1, if (m1 < P1) "no wrap" else "WRAPS",
              m2, P2, if (m2 < P2) "no wrap" else "WRAPS"))
}
cat("  so no two labels of 3 or fewer printable ASCII characters can\n")
cat("  collide: 95^1 + 95^2 + 95^3 = ",
    95 + 95^2 + 95^3, " labels, all distinct.\n", sep = "")

# and check it rather than only argue it
lab3 <- unlist(lapply(1:3, function(n) {
  g <- expand.grid(rep(list(32:126), n))
  do.call(paste0, lapply(g, intToUtf8, multiple = TRUE))
}), use.names = FALSE)
c3 <- pkg_code(lab3)
cat("  exhaustive check over all ", length(lab3),
    " of them: duplicates = ", sum(duplicated(c3)), "\n", sep = "")

# --------------------------------------------------------------------
# 2. The hunt, over identifiers of the shape a study actually assigns.
N <- 1e8
cat("\n-- 2. \"S\" plus eight digits, ", format(N, big.mark = ","),
    " labels\n", sep = "")
cat("   expected colliding pairs, birthday: ",
    sprintf("%.2f", N * (N - 1) / 2 / (P1 * P2)), "\n", sep = "")

codes <- numeric(N)
chunk <- 2e6
t0 <- Sys.time()
for (off in seq(0, N - 1, by = chunk)) {
  i <- off + seq_len(min(chunk, N - off)) - 1
  h1 <- rep(83, length(i))            # "S"
  h2 <- (1 * 137 + 83) %% P2
  h2 <- rep(h2, length(i))
  p <- i
  for (j in 7:0) {                    # the eight digits, most significant first
    d <- (p %/% 10^j) %% 10 + 48
    h1 <- (h1 * 131 + d) %% P1
    h2 <- (h2 * 137 + d) %% P2
  }
  codes[i + 1] <- h1 * 16777216 + h2
}
cat("   hashed in ", sprintf("%.1f", as.numeric(difftime(Sys.time(), t0,
                                                         units = "secs"))),
    " s\n", sep = "")

# spot check three of them against the package
sp <- c(1L, 12345678L, N)
lab_of <- function(i) sprintf("S%08d", i - 1)
cat("   spot check against the package: ",
    identical(codes[sp], pkg_code(lab_of(sp))), "\n", sep = "")

dup <- which(duplicated(codes))
cat("   colliding labels found: ", length(dup), "\n", sep = "")
if (length(dup)) {
  for (j in utils::head(dup, 6L)) {
    i <- match(codes[j], codes)
    a <- lab_of(i)
    b <- lab_of(j)
    cat(sprintf("   %s and %s both code %.0f  (package: %.0f, %.0f)\n",
                a, b, codes[j], pkg_code(a), pkg_code(b)))
  }
  saveRDS(data.frame(a = lab_of(match(codes[dup], codes)),
                     b = lab_of(dup), code = codes[dup]),
          "C:/Users/adf44/source/r/frmtmb-wt-ndt/dev/rev-ndt-collide.rds")
}
