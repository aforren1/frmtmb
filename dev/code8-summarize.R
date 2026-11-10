# Lane code8: counts for the findings, generated from the batteries' TSVs.
#   Rscript dev/code8-summarize.R
rd <- function(f) {
  p <- file.path("dev/code8-out", f)
  if (file.exists(p)) read.delim(p, quote = "", comment.char = "") else NULL
}
fmt <- function(x) if (!length(x) || all(is.na(x))) "none" else
  paste(format(range(x, na.rm = TRUE), digits = 3), collapse = " to ")
bat <- function(lab) {
  X <- rd(sprintf("battery-%s.tsv", lab))
  if (is.null(X)) return(invisible(NULL))
  c8 <- X[grepl("false", X$msg), ]
  cat(sprintf("battery %s: %d fits, %d at code 8; gain %s; dll %s;",
              lab, nrow(X), nrow(c8), fmt(c8$gain), fmt(c8$dll)),
      sprintf("se all finite %d; warned (conv_code) %d; message %d;",
              sum(c8$se_fin), sum(c8$w_code > 0), sum(c8$m_atmax > 0)),
      sprintf("with another warning %d; gain > 1e-6: %d, > 5e-6: %d\n",
              sum(c8$w_other > 0), sum(c8$gain > 1e-6), sum(c8$gain > 5e-6)))
  invisible(c8)
}
cat("== correct batteries (stops at the maximum)\n")
a <- rbind(bat("before-base"), bat("before-ywx"))
b <- rbind(bat("after-base"), bat("after-ywx"))
b2 <- rbind(bat("after-base2"), bat("after-ywx2"))
if (!is.null(b2)) cat("after (fix round), all:", nrow(b2), "code-8 fits, gain", fmt(b2$gain),
                      ", message", sum(b2$m_atmax > 0), ", warned",
                      sum(b2$w_code > 0), "
")
if (!is.null(a)) cat("before, all:", nrow(a), "code-8 fits, gain", fmt(a$gain),
                     "\n")
if (!is.null(b)) cat("after, all:", nrow(b), "code-8 fits, gain", fmt(b$gain),
                     ", message", sum(b$m_atmax > 0), ", warned",
                     sum(b$w_code > 0), "\n")
cat("== short stops\n")
for (lab in c("before", "after")) {
  X <- rd(sprintf("short-%s.tsv", lab))
  if (is.null(X)) next
  c8 <- X[grepl("false", X$msg), ]
  sh <- c8[!is.na(c8$dll) & c8$dll < -1e-6, ]
  cat(sprintf(paste("short %s: %d fits, %d at code 8, %d of them short by",
                    "more than 1e-6 (dll %s), gain on those %s;"),
              lab, nrow(X), nrow(c8), nrow(sh), fmt(sh$dll), fmt(sh$gain)),
      sprintf("warned %d of %d short, message on %d short;",
              sum(sh$w_code > 0), nrow(sh), sum(sh$m_atmax > 0)),
      sprintf(" code-8 not short: %d (gain %s, message %d)\n",
              nrow(c8) - nrow(sh), fmt(c8$gain[!(c8$case %in% sh$case)]),
              sum(c8$m_atmax[!(c8$case %in% sh$case)] > 0)))
}
cat("== ridge (Part C)\n")
for (lab in c("before", "after2")) {
  X <- rd(if (lab == "before") "ridge-before.tsv" else
    sprintf("ridge-exp-%s.tsv", lab))
  if (is.null(X)) next
  cat("ridge", lab, ":", nrow(X), "draws; objective",
      fmt(X$obj), "; se(a) finite", sum(is.finite(X$se_a)), "\n")
  tb <- table(X$lost)
  for (i in seq_along(tb)) cat("  ", tb[[i]], ":", names(tb)[i], "\n")
}
