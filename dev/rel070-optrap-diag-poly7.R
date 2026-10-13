# Release 0.70.0: the degree-7 rows of dev/optrap-p1-polystudy.R (same
# data, same seeds) with every warning's text, so a change in the warned
# count between lane optrap's build 27 and the merged build can be
# traced to the lane whose code gives or drops the warning.
#   Rscript dev/rel070-optrap-diag-poly7.R <lib or none> <label>
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "none") a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "frmtmb",
    as.character(packageVersion("frmtmb")), "label", a[2], "\n")
deg <- 7
for (s in 1:20) {
  set.seed(1000 * deg + s)
  d <- data.frame(x = runif(200, 1, 2))
  d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + ")))
  w <- character(); m <- character()
  f <- withCallingHandlers(frm(fo, data = d),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    },
    message = function(x) {
      m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
    })
  l <- lm(fo, data = d)
  gap <- as.numeric(logLik(f)) - as.numeric(logLik(l))
  se <- suppressWarnings(fixef(f))[, "Est.Error"]
  ref <- sqrt(diag(vcov(l)) * (200 - 8) / 200)
  lost <- frmtmb:::sdr_of(f)$se_lost
  cat(sprintf(paste0("seed %2d gap %10.3g code %d | lost %s | max |SE/lm-1|",
                     " %.3g | warn %d: %s | msg %d: %s\n"),
              s, gap, f$opt$convergence,
              if (length(lost)) paste(names(lost), lost, sep = ":",
                                      collapse = ",") else "none",
              max(abs(se / ref - 1), na.rm = TRUE), length(w),
              substr(gsub("[\t\n]", " ", paste(w, collapse = " || ")), 1,
                     160),
              length(m),
              substr(gsub("[\t\n]", " ", paste(m, collapse = " || ")), 1,
                     160)))
}
