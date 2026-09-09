# rev-ndt round 2, attack 1 continued: a CONSTRUCTED collision.
#
# The brute-force hunt over 10^8 identifiers of the shape "S" plus eight
# digits found none where the birthday count predicted nine, and the
# reason is not luck. A collision is a nonzero integer vector of
# per-position code-point DIFFERENCES lying in the lattice
#
#   { d : sum d_k 131^(n-k) = 0 mod P1 and sum d_k 137^(n-k) = 0 mod P2 }
#
# whose determinant is P1 * P2 = 5.6e14. A family whose labels have L
# varying positions over an alphabet spanning s code points offers a
# difference box of volume (2s+1)^L, and the expected number of nonzero
# lattice points in it is that volume over 5.6e14. Eight digits give
# 19^8 / 5.6e14 = 3e-05, so that family has no collisions to find.
#
# To find one the box has to be big. Four positions over a span of 6401
# give 6401^4 / 5.6e14 = 2.98 expected, and four positions split two and
# two is a meet in the middle rather than a search.
#
# No seed: the search is exhaustive over the stated difference box.

.libPaths(c("C:/Users/adf44/source/r/rev-ndt-lib2",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.eam))
pkg_code <- get("ddm_label_code", envir = asNamespace("frmtmb.eam"))

P1 <- 33554393
P2 <- 16777213
S <- 3200                      # |difference| at each of four positions
cat("difference box (2S+1)^4 / (P1*P2) = ",
    sprintf("%.2f", (2 * S + 1)^4 / (P1 * P2)), " expected solutions\n",
    sep = "")

A <- c(131^3, 131^2, 131, 1) %% P1
B <- c(137^3, 137^2, 137, 1) %% P2

d <- -S:S
n <- length(d)
# left half: positions 1 and 2, packed as h1 * P2 + h2, both exact
d1 <- rep(d, times = n)
d2 <- rep(d, each = n)
L <- ((d1 * A[1] + d2 * A[2]) %% P1) * P2 +
  ((d1 * B[1] + d2 * B[2]) %% P2)
cat("left keys: ", length(L), "\n", sep = "")

d3 <- d1
d4 <- d2
R <- ((-(d3 * A[3] + d4 * A[4])) %% P1) * P2 +
  ((-(d3 * B[3] + d4 * B[4])) %% P2)
cat("right keys: ", length(R), "\n", sep = "")

j <- match(L, R)
hit <- which(!is.na(j))
cat("meet-in-the-middle hits: ", length(hit), "\n", sep = "")

report <- function(i) {
  k <- j[i]
  dd <- c(d1[i], d2[i], d3[k], d4[k])
  if (all(dd == 0)) return(FALSE)
  base <- 12000L
  cb <- rep(base, 4L)
  ca <- base + dd
  if (any(ca < 33L) || any(ca > 1114111L) ||
      any(ca >= 55296L & ca <= 57343L)) return(FALSE)
  a <- intToUtf8(ca)
  b <- intToUtf8(cb)
  if (identical(a, b)) return(FALSE)
  ka <- pkg_code(a)
  kb <- pkg_code(b)
  cat("\n  difference vector: ", paste(dd, collapse = " "), "\n", sep = "")
  cat("  label A code points: ", paste(ca, collapse = " "), "\n", sep = "")
  cat("  label B code points: ", paste(cb, collapse = " "), "\n", sep = "")
  cat("  distinct labels: ", !identical(a, b), "\n", sep = "")
  cat("  ddm_label_code(A) = ", sprintf("%.0f", ka), "\n", sep = "")
  cat("  ddm_label_code(B) = ", sprintf("%.0f", kb), "\n", sep = "")
  cat("  COLLIDE: ", identical(ka, kb), "\n", sep = "")
  if (identical(ka, kb)) {
    saveRDS(list(a = a, b = b, ca = ca, cb = cb, code = ka),
            paste0("C:/Users/adf44/source/r/frmtmb-wt-ndt/dev/",
                   "rev-ndt-collide2.rds"))
    return(TRUE)
  }
  FALSE
}

found <- FALSE
for (i in hit) {
  if (report(i)) { found <- TRUE; break }
}
cat("\ncollision constructed: ", found, "\n", sep = "")
