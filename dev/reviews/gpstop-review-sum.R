# Review of lane gpstop: summary of the sweep logs (ROW lines), lane
# against norelin per row. Args: log files.
fs <- commandArgs(TRUE)
rd <- function(f) {
  l <- grep("^ROW ", readLines(f), value = TRUE)
  l <- l[!grepl(" ERROR ", l)]
  m <- regmatches(l, regexec(paste0(
    "^ROW (\\S+) (\\S+) ll (\\S+) conv (\\S+) nwarn (\\S+) wall (\\S+) ",
    "relin (\\S+) hit (\\S+) runs (\\S+) gmax (\\S+) evmin (\\S+) left (\\S+)",
    " theta (\\S+) \\| ?(.*)$"), l))
  do.call(rbind, lapply(m, function(v) data.frame(
    tag = v[2], arm = v[3], ll = as.numeric(v[4]), conv = as.integer(v[5]),
    nwarn = as.integer(v[6]), wall = as.numeric(v[7]),
    relin = as.integer(v[8]), hit = as.integer(v[9]),
    gmax = as.numeric(v[11]), evmin = as.numeric(v[12]), left = v[13],
    theta = v[14], warn = v[15])))
}
d <- do.call(rbind, lapply(fs, rd))
la <- d[d$arm == "lane", ]
no <- d[d$arm == "norelin", ]
m <- merge(la, no, by = "tag", suffixes = c(".l", ".n"))
m$dll <- m$ll.l - m$ll.n
cat("rows", nrow(m), "| relin ran", sum(m$relin.l > 0), "| relin read a",
    "new scale", sum(m$hit.l > 0), "| left closure", sum(m$left.l == "TRUE"),
    "\n")
cat("up >1e-6", sum(m$dll > 1e-6), "up >1e-3", sum(m$dll > 1e-3),
    "down < -1e-6", sum(m$dll < -1e-6), "max up", max(m$dll), "\n")
cat("warn or code changes", sum(m$nwarn.l != m$nwarn.n | m$conv.l != m$conv.n),
    "\n")
th_ext <- function(s) {
  v <- as.numeric(strsplit(s, ",")[[1]])
  any(abs(v) > 8)
}
m$deg.l <- m$evmin.l < 1e-4 | m$gmax.l > 1e-2 |
  vapply(m$theta.l, th_ext, NA)
m$deg.n <- m$evmin.n < 1e-4 | m$gmax.n > 1e-2 |
  vapply(m$theta.n, th_ext, NA)
cat("suspect point (evmin<1e-4, gmax>1e-2 or |theta|>8): lane", sum(m$deg.l),
    "norelin", sum(m$deg.n), "lane-only", sum(m$deg.l & !m$deg.n), "\n")
cat(sprintf("wall lane %.1f norelin %.1f ratio %.3f | rows relin hit: ratio %.3f | rows no hit: ratio %.3f\n",
            sum(m$wall.l), sum(m$wall.n), sum(m$wall.l) / sum(m$wall.n),
            sum(m$wall.l[m$hit.l > 0]) / sum(m$wall.n[m$hit.l > 0]),
            sum(m$wall.l[m$hit.l == 0]) / sum(m$wall.n[m$hit.l == 0])))
show <- m[abs(m$dll) > 1e-6 | m$deg.l | m$deg.n |
            m$nwarn.l != m$nwarn.n | m$conv.l != m$conv.n, ]
for (i in seq_len(nrow(show))) {
  s <- show[i, ]
  cat(sprintf(paste0("SUMROW %s dll %+.4f | lane %.4f c%d w%d g %.1e ev %.1e ",
                     "th %s %s | norelin %.4f c%d w%d g %.1e ev %.1e th %s %s\n"),
              s$tag, s$dll, s$ll.l, s$conv.l, s$nwarn.l, s$gmax.l, s$evmin.l,
              s$theta.l, s$warn.l, s$ll.n, s$conv.n, s$nwarn.n, s$gmax.n,
              s$evmin.n, s$theta.n, s$warn.n))
}
