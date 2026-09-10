# Lane `latent`, punch round 2: the shrink grid, in sample and out,
# NEVER pooled.
#
#   Rscript dev/latent-2p4-weight-summarize.R
#
# In sample  : the lane's own 200 seeds, 20260910 to 20261109.
#              w = 0, 0.1, 0.25, 0.5, 0.75, 0.9 from
#              dev/latent-2p4-shrink.tsv; w = 0.95, 0.99, 1 from
#              dev/latent-2p4-shrink2.tsv; the 0.2.2 score cut from
#              dev/latent-2p4-lca.tsv.
# Out of sample: 200 fresh seeds, 20270401 to 20270600, from
#              dev/latent-2p4-oos.tsv.
#
# THE REFERENCE IS ADJUDICATED. `poLCA(nrep = 10)` can miss too, so a
# seed's reference is the LARGEST log-likelihood anything reached on it:
# poLCA's best of ten, or any start in the grid. A start LOSES a seed
# when it is more than 1e-6 relative below that. Taking poLCA on trust
# would flatter every start that agrees with it and punish one that
# beats it.

source("dev/latent-env.R")

lost_of <- function(lls, ref, ll_names) {
  tol <- 1e-6 * abs(ref)
  vapply(ll_names, function(nm) sum(lls[[nm]] < ref - tol, na.rm = TRUE),
         integer(1))
}

## ---------------------------------------------------------- in sample
a <- utils::read.delim("dev/latent-2p4-shrink.tsv", check.names = FALSE)
b <- utils::read.delim("dev/latent-2p4-shrink2.tsv", check.names = FALSE)
z <- utils::read.delim("dev/latent-2p4-lca.tsv", check.names = FALSE)
stopifnot(identical(a$seed, b$seed), identical(a$seed, z$seed))
ins <- data.frame(seed = a$seed, polca = a$ll_polca, cut = z$ll_frm)
for (w in c("0", "01", "025", "05", "075", "09")) {
  ins[[paste0("w", w)]] <- a[[paste0("gap_w", w)]] + a$ll_polca
}
ins[["w095"]] <- b[["ll_w095"]]
ins[["w099"]] <- b[["ll_w099"]]
ins[["w1"]] <- b[["ll_w100"]]
wn <- c("cut", "w0", "w01", "w025", "w05", "w075", "w09", "w095",
        "w099", "w1")
ref_in <- do.call(pmax, c(list(ins$polca), ins[wn], list(na.rm = TRUE)))
lost_in <- lost_of(ins, ref_in, wn)
beat_in <- vapply(wn, function(nm)
  sum(ins[[nm]] > ins$polca + 1e-6 * abs(ins$polca), na.rm = TRUE),
  integer(1))

## ------------------------------------------------------ out of sample
o <- utils::read.delim("dev/latent-2p4-oos.tsv", check.names = FALSE)
oo <- data.frame(seed = o$seed, polca = o$ll_polca)
map <- c(w0 = "ll_w000", w025 = "ll_w025", w05 = "ll_w050",
         w075 = "ll_w075", w09 = "ll_w090", w095 = "ll_w095",
         w099 = "ll_w099", w1 = "ll_w100")
for (nm in names(map)) oo[[nm]] <- o[[map[[nm]]]]
wn2 <- names(map)
ref_out <- do.call(pmax, c(list(oo$polca), oo[wn2], list(na.rm = TRUE)))
lost_out <- lost_of(oo, ref_out, wn2)
beat_out <- vapply(wn2, function(nm)
  sum(oo[[nm]] > oo$polca + 1e-6 * abs(oo$polca), na.rm = TRUE),
  integer(1))

lab <- c(cut = "the 0.2.2 score cut", w0 = "0", w01 = "0.1",
         w025 = "0.25", w05 = "0.5", w075 = "0.75", w09 = "0.9",
         w095 = "0.95", w099 = "0.99", w1 = "1")
cat("in sample : ", nrow(ins), " seeds, ", min(ins$seed), " to ",
    max(ins$seed), "\n", sep = "")
cat("out of sample: ", nrow(oo), " seeds, ", min(oo$seed), " to ",
    max(oo$seed), "\n\n", sep = "")
cat(sprintf("  %-20s %12s %14s\n", "w", "in-sample", "out-of-sample"))
for (nm in wn) {
  io <- if (nm %in% wn2) sprintf("%d / %d", lost_out[[nm]], nrow(oo)) else
    "not run"
  cat(sprintf("  %-20s %7d / %d %14s\n", lab[[nm]], lost_in[[nm]],
              nrow(ins), io))
}
cat("\n  seeds where poLCA(nrep = 10) was itself beaten, in sample:",
    sum(beat_in), " out of sample:", sum(beat_out), "\n")

cat("\n== the seeds each candidate loses ==\n")
for (nm in c("cut", "w05", "w075", "w09", "w095")) {
  ii <- ins$seed[which(ins[[nm]] < ref_in - 1e-6 * abs(ref_in))]
  oi <- if (nm %in% wn2) {
    oo$seed[which(oo[[nm]] < ref_out - 1e-6 * abs(ref_out))]
  } else integer(0)
  cat("  ", lab[[nm]], ": in ", paste(ii, collapse = " "),
      " | out ", paste(oi, collapse = " "), "\n", sep = "")
}

cat("\n== how deep the losses are, in log-likelihood units ==\n")
for (nm in c("cut", "w05", "w075", "w09", "w095")) {
  d <- c(ref_in - ins[[nm]],
         if (nm %in% wn2) ref_out - oo[[nm]] else numeric(0))
  d <- d[is.finite(d) & d > 1e-6 * abs(c(ref_in, if (nm %in% wn2)
    ref_out else numeric(0)))]
  cat("  ", lab[[nm]], ": ",
      if (length(d)) formatC(sort(d), digits = 6, format = "g") else
        "none", "\n", sep = " ")
}

cat("\n== the symmetry axis, for the upper edge ==\n")
cat("  w = 1 loses", lost_in[["w1"]], "of", nrow(ins), "in sample and",
    lost_out[["w1"]], "of", nrow(oo), "out of sample\n")
cat("  its worst gap:",
    formatC(max(c(ref_in - ins$w1, ref_out - oo$w1), na.rm = TRUE),
            digits = 6, format = "g"), "units\n")
