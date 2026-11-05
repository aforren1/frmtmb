# Reads the CALIB lines of dev/gpgap-calib.R and scores candidate rules for the
# nugget warning against the reference.
#   Rscript dev/gpgap-calib-sum.R [print]
L <- unlist(lapply(Sys.glob("dev/gpgap-out/calib/c*.txt"), readLines, warn = FALSE))
L <- L[startsWith(L, "CALIB")]
# a rerun of a row replaces it; the file sort puts reruns last
L <- L[!duplicated(vapply(strsplit(L, " | ", fixed = TRUE), `[`, "", 1),
                    fromLast = TRUE)]
bad <- L[grepl("ERROR|REF FAILED", L)]
L <- L[!grepl("ERROR|REF FAILED", L)]
num <- function(pat) suppressWarnings(as.numeric(sub(pat, "\\1", L)))
tb <- data.frame(
  tag = sub("^CALIB (.*?) \\| m .*", "\\1", L, perl = TRUE),
  m = num(".*\\| m (\\d+) rows.*"),
  rows = num(".* rows (\\d+) conv.*"),
  conv = num(".* rows \\d+ conv (\\d+) warned.*"),
  warned = grepl("warned TRUE", L),
  rmed = num(".*rmed (\\S+) rmean.*"),
  rmean = num(".*rmean (\\S+) rq90.*"),
  rq90 = num(".*rq90 (\\S+) rmax.*"),
  rmax = num(".*rmax (\\S+) scan10.*"),
  scan = num(".*scan10 (\\S+) \\|.*"),
  sd_f = num(".*sd_f (\\S+) sd_r.*"),
  sd_r = num(".*sd_r (\\S+) \\|.*"),
  gap = num(".* gap (\\S+) \\|.*"),
  dz = num(".*dzmax (\\S+) dzsd.*"),
  dzsd = num(".*dzsd (\\S+) \\|.*"),
  refconv = num(".*refconv (\\d+) refgr.*"),
  refgr = num(".*refgr (\\S+) \\|.*"),
  stringsAsFactors = FALSE)
tb$fam <- sub(" .*", "", tb$tag)
# Effect: the reference's Wald reading of how far frm() sits. ll gap first:
# 0.5 nats is a quarter of the one-parameter 95 percent boundary.
tb$moved <- !is.na(tb$dz) & tb$dz >= 1
tb$big <- abs(tb$gap) >= 0.5
# far above the reference: it stopped short (a nugget can pay a few tenths of a
# nat in a model with known noise, dev/gpgap-binom1e5.R, but not this)
tb$unread <- tb$gap < -2
rules <- list(
  old = tb$warned,
  r_med = tb$rmed >= sqrt(2 / tb$m),
  r_mean = tb$rmean >= sqrt(2 / tb$m),
  r_q90 = tb$rq90 >= sqrt(2 / tb$m),
  scan = tb$scan >= 1,
  fix_0.1 = tb$rmed >= 0.1,
  fix_0.25 = tb$rmed >= 0.25
)
score <- function(r, ok = !tb$unread) {
  mv <- tb$moved | tb$big
  sprintf("fires %3d | of %3d moved: hit %3d miss %3d | of %3d unmoved: false %3d",
          sum(r & ok, na.rm = TRUE), sum(mv & ok), sum(r & mv & ok, na.rm = TRUE),
          sum(!r & mv & ok, na.rm = TRUE), sum(!mv & ok),
          sum(r & !mv & ok, na.rm = TRUE))
}
cat("rows", nrow(tb), "unreadable", sum(tb$unread), "errors", length(bad), "\n")
cat(substr(bad, 1, 150), sep = "\n")
for (nm in names(rules)) cat(sprintf("%-9s %s\n", nm, score(rules[[nm]])))
if (length(commandArgs(trailingOnly = TRUE))) {
  options(width = 180)
  tb$r <- rules$r_q90
  o <- tb[, c("tag", "m", "conv", "warned", "rmed", "rq90", "scan", "sd_f", "sd_r",
              "gap", "dz", "dzsd", "refconv", "refgr", "moved", "big", "r")]
  o[] <- lapply(o, function(x) if (is.numeric(x)) signif(x, 3) else x)
  print(o[order(o$tag), ], row.names = FALSE)
}
