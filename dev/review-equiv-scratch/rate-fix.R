# Reviewer: the headline under the review's verdict corrections, and
# sensitivities. Same estimator as dev/equiv-summary.R (stratified, reach
# shares, normal interval with p(1-p)/n per stratum).
r <- read.delim("dev/equiv-data/results.tsv", stringsAsFactors = FALSE)
cp <- read.delim("dev/equiv-frame-out/compat.tsv")
cp <- cp[cp$in_core_only, ]
w <- read.delim("dev/equiv-data/weights.tsv")
pairw <- function(wt) wt$weight[match(cp$feature_a, wt$feature)] *
  wt$weight[match(cp$feature_b, wt$feature)]
cp$w <- pairw(w)
est <- function(agree, wnew = NULL) {
  W <- tapply(if (is.null(wnew)) cp$w else wnew, cp$status, sum)
  W <- W / sum(W)
  p <- v <- 0
  for (s in names(W)) {
    k <- r$status == s
    y <- rep(as.numeric(agree[k]), r$draws[k])
    rr <- if (is.null(wnew)) rep(1, sum(k)) else {
      key <- paste(cp$feature_a, cp$feature_b)
      m <- match(paste(r$feature_a, r$feature_b)[k], key)
      wnew[m] / cp$w[m]
    }
    rr <- rep(rr, r$draws[k])
    ph <- sum(rr * y) / sum(rr)
    neff <- sum(rr)^2 / sum(rr^2)
    p <- p + W[[s]] * ph
    v <- v + W[[s]]^2 * ph * (1 - ph) / neff
  }
  sprintf("%.1f [%.1f, %.1f]", 100 * p, 100 * (p - 1.96 * sqrt(v)),
          100 * (p + 1.96 * sqrt(v)))
}
show <- function(lab, a) cat(sprintf("%-64s %s\n", lab, a))
show("A lane, as reported", est(r$agree))
B <- r$agree
B[r$id %in% c("p158", "p167", "p029")] <- FALSE
show("B p158, p167 -> divergence by decision; p029 -> cannot measure",
     est(B))
C <- B
C[r$id %in% c("p042", "p047", "p143", "p136", "p044")] <- FALSE
show("C = B, and L4-posterior-only penalized pairs and REML uncounted",
     est(C))
w2 <- w
w2$weight[w2$feature == "bar_crossing"] <- 1
show("D = B, bar_crossing reach weight 3 -> 1", est(B, pairw(w2)))
show("E = C, bar_crossing reach weight 3 -> 1", est(C, pairw(w2)))
# measured only, under B
M <- r$verdict != "cannot measure" & r$id != "p029"
cat("\nmeasured-only under B:\n")
Wm <- tapply(cp$w, cp$status, sum); Wm <- Wm / sum(Wm)
p <- 0
for (s in names(Wm)) {
  k <- r$status == s & M
  p <- p + Wm[[s]] * sum(r$draws[k] * B[k]) / sum(r$draws[k])
}
cat(sprintf("  %.1f\n", 100 * p))
cat("\nrefused stratum under B:", sum(r$draws[r$status == "refused" & B]),
    "of", sum(r$draws[r$status == "refused"]), "\n")
nd <- r$nodefect
nd[r$id %in% c("p029", "p022", "p023", "p094", "p109")] <- FALSE
show("no defect: lane", est(r$nodefect))
show("no defect: p029 and the four exact gp() pairs -> cannot measure",
     est(nd))
