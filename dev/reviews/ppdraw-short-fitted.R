# Short check: fitted() standard errors on 3269 (zib, zi at its end).
# Merged build (lane + boundfam) against 0.81.0 and against the limit
# model (binomial), and the row gradient along zi's direction.
# Usage: Rscript ppdraw-short-fitted.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(7)
d <- data.frame(x = stats::rnorm(200), N = sample(5:15, 200, TRUE))
d$y <- stats::rbinom(200, d$N, stats::plogis(0.3 + 0.5 * d$x))
fz <- suppressMessages(frm(bf(y | trials(N) ~ x),
                           family = zero_inflated_binomial(), data = d))
fb <- frm(bf(y | trials(N) ~ x), family = binomial(), data = d)
w <- character(0)
fe <- withCallingHandlers(fitted(fz), warning = function(x) {
  w <<- c(w, conditionMessage(x))
  invokeRestart("muffleWarning")
})
cat("warnings:", if (length(w)) substr(w, 1, 140) else "(none)", "\n")
eb <- fitted(fb)
cat("zib fitted Est.Error: NaN rows", sum(is.nan(fe[, "Est.Error"])),
    " first rows", format(fe[1:4, "Est.Error"], digits = 4), "\n")
cat("binomial (limit model) Est.Error first rows",
    format(eb[1:4, "Est.Error"], digits = 4), "\n")
cat("max |zib - binomial| Est.Error over finite rows:",
    format(max(abs(fe[, "Est.Error"] - eb[, "Est.Error"]), na.rm = TRUE),
           digits = 4), "\n")
cat("Estimate max |zib - binomial|:",
    format(max(abs(fe[, "Estimate"] - eb[, "Estimate"])), digits = 4), "\n")
# the rows' derivative along the zi intercept, by central difference
map <- frmtmb:::outer_par_map(fz)
v0 <- frmtmb:::fit_outer_vector(fz, map)
j <- grep("^zi", map$names)
h <- 1e-3
up <- v0; up[j] <- up[j] + h
dn <- v0; dn[j] <- dn[j] - h
gu <- fitted(frmtmb:::fit_set_outer(fz, up, map))[, "Estimate"]
gd <- fitted(frmtmb:::fit_set_outer(fz, dn, map))[, "Estimate"]
gz <- (gu - gd) / (2 * h)
cat("d fitted / d zi_(Intercept): max |.|", format(max(abs(gz)), digits = 4),
    "\n")
sdr <- frmtmb:::sdr_of(fz)
cat("se_lost:", paste(names(sdr$se_lost), sdr$se_lost), "\n")
if (!is.null(sdr$se_null)) {
  cat("se_null direction:\n")
  print(signif(sdr$se_null, 4))
}
