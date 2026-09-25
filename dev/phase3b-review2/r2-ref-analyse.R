# Reviewer 2, item 2: the lane against my Rmpfr reference (r2-ref.R).
# Error = absolute error of the LOG (the relative error of the quantity).
# Rows whose reference is below log(1e-300) are reported separately for
# F: frmtmb floors left-censored probabilities at 1e-300.
x <- readRDS("dev/phase3b-review2/ref-points.rds")
options(width = 170, digits = 4)
x$va <- abs(x$v * x$a)
x$eS <- abs(x$g_lS - x$lS); x$eF <- abs(x$g_lF - x$lF)
x$eFl <- abs(x$g_lFl - x$lFl); x$eFu <- abs(x$g_lFu - x$lFu)
cat("points", nrow(x), "; non-finite lane values:",
    sum(!is.finite(x$g_lS) | !is.finite(x$g_lF) | !is.finite(x$g_lFl) | !is.finite(x$g_lFu)), "\n")
cat("series agreement where both ran (u <= 3): S max", max(x$agreeS, na.rm = TRUE),
    "; F_lower where lFl > -600: max", max(x$agreeFl[x$lFl > -600], na.rm = TRUE), "\n")
band <- cut(x$va, c(-1, 1, 5, 12, 24, 48, 72, 120, 400))
rel <- function(e, ref) ifelse(abs(ref) > 1, e / abs(ref), e)
cat("\nabsolute log error by |v| a band (all rows; F columns only where ref > log 1e-300):\n")
ok <- function(r) r > log(1e-300)
tb <- data.frame(
  S = tapply(x$eS, band, max),
  F = tapply(ifelse(ok(x$lF), x$eF, NA), band, max, na.rm = TRUE),
  Fl = tapply(ifelse(ok(x$lFl), x$eFl, NA), band, max, na.rm = TRUE),
  Fu = tapply(ifelse(ok(x$lFu), x$eFu, NA), band, max, na.rm = TRUE),
  S_relLog = tapply(rel(x$eS, x$lS), band, max),
  rows = as.vector(table(band)))
print(tb)
cat("\nworst 12 on S:\n")
print(head(x[order(-x$eS), c("u", "w", "a", "v", "va", "lS", "g_lS", "eS")], 12))
cat("\nworst 12 on F_lower/F_upper (ref > log 1e-300):\n")
y <- rbind(data.frame(x[ok(x$lFl), c("u", "w", "a", "v", "va")], ref = x$lFl[ok(x$lFl)],
                      got = x$g_lFl[ok(x$lFl)], err = x$eFl[ok(x$lFl)], b = "lower"),
           data.frame(x[ok(x$lFu), c("u", "w", "a", "v", "va")], ref = x$lFu[ok(x$lFu)],
                      got = x$g_lFu[ok(x$lFu)], err = x$eFu[ok(x$lFu)], b = "upper"))
print(head(y[order(-y$err), ], 12))
cat("\nworst 8 on F (both):\n")
print(head(x[ok(x$lF), ][order(-x$eF[ok(x$lF)]), c("u", "w", "a", "v", "va", "lF", "g_lF", "eF")], 8))
cat("\nrows with the reference below log 1e-300 (F lower/upper): lane error there, max",
    max(c(x$eFl[!ok(x$lFl)], x$eFu[!ok(x$lFu)])), "(relative:",
    max(c(x$eFl[!ok(x$lFl)] / abs(x$lFl[!ok(x$lFl)]), x$eFu[!ok(x$lFu)] / abs(x$lFu[!ok(x$lFu)]))), ")\n")
cat("\nTaylor branch rows (|v a| < 0.1): max error S", max(x$eS[x$va < 0.1]),
    " Fl", max(x$eFl[x$va < 0.1 & ok(x$lFl)]), " Fu", max(x$eFu[x$va < 0.1 & ok(x$lFu)]), "\n")
cat("u <= 1e-4 rows: max error S", max(x$eS[x$u <= 1e-4]), " F", max(x$eF[x$u <= 1e-4 & ok(x$lF)]), "\n")
cat("u >= 40 rows: max error S", max(x$eS[x$u >= 40]), "; smallest ref lS", min(x$lS), "\n")
cat("a = 0.1 rows: max S", max(x$eS[x$a == 0.1]), "; a = 10: ", max(x$eS[x$a == 10]), "\n")
cat("w <= 0.001 or >= 0.999: max S", max(x$eS[x$w <= 0.001 | x$w >= 0.999]), "\n")
