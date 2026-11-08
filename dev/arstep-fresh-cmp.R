# Lane arstep, fix round: the reviewer's 168 fresh rows on the fixed build.
# Reads the reviewer's logs (dev/review-arstep-out/fresh-<part>-base.txt and
# -lane.txt, whose FRR lines are its multi-start references) and this
# lane's rerun (dev/arstep-out/fresh-<part>-fix.txt, review-arstep-fresh.R
# lane <part> 0). Reference per row: the best of the base, the first lane
# build, its references and the fixed build. Short: more than 1e-3 under it.
# Worse: the fixed build under the base by more than 1e-3. Args: none
rd <- function(f) {
  L <- readLines(f)
  fr <- grep("^FR ", L, value = TRUE)
  fr <- fr[!grepl(" ERROR ", fr)]
  sp <- strsplit(fr, " +")
  x <- data.frame(i = as.integer(vapply(sp, `[`, "", 3)),
                  ll = as.numeric(vapply(sp, function(s) s[which(s == "ll") + 1], "")),
                  conv = as.integer(vapply(sp, function(s) s[which(s == "conv") + 1], "")),
                  w = as.integer(vapply(sp, function(s) s[which(s == "w") + 1], "")),
                  wall = as.numeric(vapply(sp, function(s) s[which(s == "wall") + 1], "")))
  rr <- grep("^FRR ", L, value = TRUE)
  ref <- if (length(rr)) {
    s2 <- strsplit(rr, " +")
    setNames(as.numeric(vapply(s2, function(s) s[which(s == "best") + 1], "")), vapply(s2, `[`, "", 3))
  } else numeric(0)
  list(x = x, ref = ref, lib = grep("^lib:", L, value = TRUE), err = sum(grepl(" ERROR ", L)))
}
rv <- "C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/review-arstep-out"
me <- "C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-out"
tot <- c(rows = 0, bshort = 0, fshort = 0, worse = 0, better = 0)
for (part in c("G", "N", "P")) {
  b <- rd(file.path(rv, sprintf("fresh-%s-base.txt", part)))
  l <- rd(file.path(rv, sprintf("fresh-%s-lane.txt", part)))
  ff <- file.path(me, sprintf("fresh-%s-fix.txt", part))
  if (!file.exists(ff)) next
  f <- rd(ff)
  cat("FC", part, "logs:", b$lib, "|", f$lib, "| errors base", b$err, "fix", f$err, "\n")
  ids <- intersect(intersect(b$x$i, l$x$i), f$x$i)
  g <- function(d, i) d$x[match(i, d$x$i), ]
  B <- g(b, ids); Lx <- g(l, ids); Fx <- g(f, ids)
  ref <- pmax(B$ll, Lx$ll, Fx$ll, l$ref[as.character(ids)], b$ref[as.character(ids)], na.rm = TRUE)
  bs <- B$ll - ref < -1e-3; fs <- Fx$ll - ref < -1e-3
  worse <- Fx$ll < B$ll - 1e-3; better <- Fx$ll > B$ll + 1e-3
  cat(sprintf("FC %s rows %d | base short %d worst %.3g | fixed short %d silent %d worst %.3g | fixed worse than base %d, better %d | wall base %.1f fixed %.1f (not side by side)\n",
              part, length(ids), sum(bs), if (any(bs)) min(B$ll - ref) else 0, sum(fs), sum(fs & Fx$conv == 0 & Fx$w == 0),
              if (any(fs)) min(Fx$ll - ref) else 0, sum(worse), sum(better), sum(B$wall), sum(Fx$wall)))
  if (any(worse)) cat("FC worse rows:", paste(sprintf("%d (%.3g)", ids[worse], (Fx$ll - B$ll)[worse]), collapse = "; "), "\n")
  dl <- Fx$ll - Lx$ll
  cat(sprintf("FC %s fixed against the first lane build: lower by >1e-3 on %d, higher on %d\n", part, sum(dl < -1e-3), sum(dl > 1e-3)))
  tot <- tot + c(length(ids), sum(bs), sum(fs), sum(worse), sum(better))
}
cat(sprintf("FCSUM rows %d base short %d fixed short %d fixed worse than base %d better %d\n",
            tot[1], tot[2], tot[3], tot[4], tot[5]))
cat("DONE\n")
