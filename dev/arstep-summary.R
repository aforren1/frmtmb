# Lane arstep: counts for the findings, generated from the logs in
# dev/arstep-out. Reference per row: the best log-likelihood any arm or any
# reference start reached on that row (the tape is the same on every arm).
# Short: more than 1e-3 nats under the reference. Silent: conv 0, no
# warning. Wall: the default fit's, summed over the rows every listed arm has.
# Args: none.
out <- "C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-out"
rd <- function(f, tag) {
  p <- file.path(out, f)
  if (!file.exists(p)) return(NULL)
  l <- grep(paste0("^", tag, " "), readLines(p), value = TRUE)
  l <- l[!grepl("ERROR", l)]
  if (!length(l)) return(NULL)
  l
}
num <- function(l, key) as.numeric(sub(paste0(".* ", key, " ([^ ]+).*"), "\\1", l))
# arms: named list of log lines; refs: extra lines whose ll and best enter
# the reference only; base: the arm the others are compared with
tab <- function(arms, refs, label, keycols, base = names(arms)[1L]) {
  key <- function(l) vapply(strsplit(l, " +"), function(v) paste(v[keycols], collapse = " "), "")
  ref <- numeric(0)
  add <- function(k, v) {
    for (i in seq_along(k)) if (is.finite(v[i])) ref[k[i]] <<- max(ref[k[i]], v[i], na.rm = TRUE)
  }
  for (l in c(arms, refs)) {
    if (!length(l)) next
    k <- key(l)
    add(k, num(l, "ll"))
    if (any(grepl(" best ", l))) add(k[grepl(" best ", l)], num(l[grepl(" best ", l)], "best"))
    if (any(grepl(" cf ", l))) add(k[grepl(" cf ", l)], num(l[grepl(" cf ", l)], "cf"))
  }
  common <- Reduce(intersect, lapply(arms, key))
  vals <- list()
  for (nm in names(arms)) {
    l <- arms[[nm]]
    k <- key(l); ll <- num(l, "ll"); g <- ll - ref[k]
    conv <- num(l, "conv"); w <- num(l, "w"); wall <- num(l, "wall")
    sh <- g < -1e-3
    cat(sprintf("SUM %-26s %-10s rows %3d short %3d silent %3d worst %9.3g | wall on %d common rows %8.1f s\n",
                label, nm, length(l), sum(sh), sum(sh & conv == 0 & w == 0), if (any(sh)) min(g) else 0,
                length(common), sum(wall[k %in% common])))
    vals[[nm]] <- setNames(ll, k)
  }
  b <- vals[[base]]
  for (nm in setdiff(names(arms), base)) {
    v <- vals[[nm]]
    cm <- intersect(names(b), names(v))
    worse <- cm[v[cm] < b[cm] - 1e-3]
    better <- cm[v[cm] > b[cm] + 1e-3]
    cat(sprintf("SUM %-26s %-10s against %s: worse %d, better %d%s\n", label, nm, base, length(worse), length(better),
                if (length(worse)) paste0(": ", paste(sprintf("%s (%.3g)", worse, v[worse] - b[worse]), collapse = "; ")) else ""))
    if (length(better)) cat(sprintf("SUMB %-26s %-10s better: %s\n", label, nm,
                                    paste(sprintf("%s (+%.3g)", better, v[better] - b[better]), collapse = "; ")))
  }
  sh_rows <- function(l) { k <- key(l); g <- num(l, "ll") - ref[k]; k[g < -1e-3] }
  for (nm in names(arms)) {
    s <- sh_rows(arms[[nm]])
    if (length(s)) {
      k <- key(arms[[nm]]); g <- num(arms[[nm]], "ll") - ref[k]
      cat(sprintf("SUMS %-26s %-10s short rows: %s\n", label, nm,
                  paste(sprintf("%s (%.3g)", k[k %in% s], g[k %in% s]), collapse = "; ")))
    }
  }
}
sub_rows <- function(l, pat) l[grepl(pat, l)]

# 150 off-gaussian rows: dev/arstep-ng.R (ngt-*: default fit only, base and
# lane run side by side; ng-base-*: the base with the 15 reference starts;
# ng-lane-none / -local: two strategies on an earlier build of the lane)
# fix round: the fixed build is the lane arm, the first build a reference
fx <- function(f, g) if (file.exists(file.path(out, f))) f else g
arms_ng <- list(base = rd("ngt-base.txt", "ANG"),
                lane = rd(fx("ngt-lane-fix.txt", "ngt-lane.txt"), "ANG"))
refs_ng <- list(rd("ngt-lane.txt", "ANG"), rd("ng-base-a.txt", "ANG"), rd("ng-base-b.txt", "ANG"),
                rd("ng-lane-none.txt", "ANG"), rd("ng-lane-local.txt", "ANG"))
if (all(lengths(arms_ng) > 0)) {
  for (grp in list(c("all", "."), c("a: bern pois binom5", " (bern|pois_low|pois_mid|binom5) "),
                   c("b: the other six", " (negbin|gamma|lognorm|student|gsigma|gar) "),
                   c("step curve", " step "), c("student", " student "), c("gar", " gar "))) {
    tab(lapply(arms_ng, sub_rows, grp[2]), lapply(refs_ng, sub_rows, grp[2]), paste0("ng ", grp[1]), 3:6)
  }
  nb <- c(refs_ng[[1]], refs_ng[[2]])
  cat(sprintf("SUM ng reference starts ran on %d of 150 base rows\n", length(nb)))
}

# 192 gaussian curves: dev/arstep-gauss.R, reference the closed form
arms_g <- list(base = rd("gauss-base.txt", "AG"),
               lane = rd(fx("gauss-lane-fix.txt", "gauss-lane.txt"), "AG"))
if (all(lengths(arms_g) > 0)) tab(arms_g, list(), "gauss 192", 3:6)

# 72 rows of dev/review-ncgp-r2-ng.R (dev/arstep-r2ng.R); its 48 are the
# eight spellings other than the two Poissons, bernoulli and 1000 trials
r2 <- function(f) {
  l <- rd(f, "RNG2")
  if (is.null(l)) return(NULL)
  # put the fields in the ANG layout the table reads
  sub(" default ll ", " ll ", sub(" [|] best ", " best ", sub(" ([0-9.]+)s$", " wall \\1", l)))
}
arms_r2 <- list(base = r2("r2ng-base.txt"),
                lane = r2(fx("r2ng-lane-fix.txt", "r2ng-lane.txt")))
refs_r2 <- list(r2("r2ng-base-ref.txt"), r2("r2ng-lane.txt"))
if (all(lengths(arms_r2) > 0)) {
  tab(arms_r2, refs_r2, "r2ng 72", c(3, 4, 6, 8))
  tab(lapply(arms_r2, sub_rows, " (negbin|gamma|lognorm|student|gsigma|gmv|glog|gar) "),
      lapply(refs_r2, sub_rows, " (negbin|gamma|lognorm|student|gsigma|gmv|glog|gar) "), "r2ng the 48", c(3, 4, 6, 8))
}
cat("DONE\n")
