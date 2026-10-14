# Release reviewer of 0.71.0: the verdict moves of the merged build from
# 0.70.0 (base) over all seeds of dev/relrev071-falarm.sh runs, outside
# negbin_large_shape (shape MLE infinite): new warnings, warnings gone,
# and fits that went silent with a worse SE.
#   Rscript dev/relrev071-falarm-moves.R <dir> [<dir> ...]
dirs <- commandArgs(trailingOnly = TRUE)
rd <- function(f) utils::read.delim(f, colClasses = "character")
get <- function(cfg, arm) do.call(rbind, lapply(dirs, function(d)
  do.call(rbind, lapply(c("nl", "lin", "poly", "mixed", "bounded"),
    function(g) rd(file.path(d, "falarm",
                             sprintf("%s-%s-%s.tsv", g, cfg, arm)))))))
for (cfg in c("ref", "ob")) {
  m <- get(cfg, "merged"); b <- get(cfg, "base")
  stopifnot(identical(m$model, b$model), identical(m$seed, b$seed))
  keep <- m$model != "negbin_large_shape"
  wm <- as.integer(m$se_warn) + as.integer(m$other_warn)
  wb <- as.integer(b$se_warn) + as.integer(b$other_warn)
  rm <- suppressWarnings(as.numeric(m$max_rel))
  rb <- suppressWarnings(as.numeric(b$max_rel))
  cat("\n==", cfg, ": fits", sum(keep), "outside negbin_large_shape\n")
  sel <- keep & (wm > wb | (m$code != "0" & b$code == "0"))
  cat("  merged warns more or has a nonzero code where base had 0:",
      sum(sel), "\n")
  if (any(sel)) print(data.frame(model = m$model[sel], seed = m$seed[sel],
    code = paste(m$code[sel], b$code[sel], sep = "/"),
    warn = paste(wm[sel], wb[sel], sep = "/"),
    lost = paste(m$lost[sel], b$lost[sel], sep = " / "),
    rel = paste(m$max_rel[sel], b$max_rel[sel], sep = " / ")),
    row.names = FALSE)
  sel <- keep & wm < wb
  cat("  merged warns less:", sum(sel), "\n")
  if (any(sel)) print(data.frame(model = m$model[sel], seed = m$seed[sel],
    code = paste(m$code[sel], b$code[sel], sep = "/"),
    warn = paste(wm[sel], wb[sel], sep = "/"),
    lost = paste(m$lost[sel], b$lost[sel], sep = " / "),
    rel = paste(m$max_rel[sel], b$max_rel[sel], sep = " / ")),
    row.names = FALSE)
  silent <- keep & wm == 0 & m$code == "0"
  sel <- silent & is.finite(rm) & is.finite(rb) & rm > 2 * rb & rm > 0.01
  cat("  silent on merged with max_rel > 1% and > 2x base:", sum(sel), "\n")
  if (any(sel)) print(data.frame(model = m$model[sel], seed = m$seed[sel],
    rel = paste(m$max_rel[sel], b$max_rel[sel], sep = " / ")),
    row.names = FALSE)
}
