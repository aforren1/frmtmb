# Punch round 1, m4: of the mo() seeds far below the profile maximum,
# how many would a cheap plateau test catch? The test is
# dev/nanse-mo-escape.R's first round: a 1 percent step toward a simplex
# vertex that lowers the objective (rounds > 0 there). Read from that
# script's TSV and the lane's SE results.
#   Rscript dev/nanse-p1-m4.R
X <- utils::read.delim("dev/nanse-log/mo-escape.tsv")
L <- utils::read.delim("dev/nanse-log/mo-lane-p1.tsv")
L <- L[L$form == "int", ]
m <- merge(X, L[, c("seed", "se_finite", "n_par", "n_warn_fit")])
big <- m$gap_fit > 0.5
cat("seeds more than 0.5 below the maximum:", sum(big), "\n")
cat("  of those silent with every SE finite on the lane:",
    sum(big & m$se_finite == m$n_par & m$n_warn_fit == 0), "\n")
cat("  of those a vertex step improves (plateau test fires):",
    sum(big & m$rounds > 0), "\n")
cat("plateau test fires on seeds within 1e-3 of the maximum:",
    sum(m$gap_fit <= 1e-3 & m$rounds > 0), "of", sum(m$gap_fit <= 1e-3),
    "\n")
cat("plateau test fires in all:", sum(m$rounds > 0), "\n")
