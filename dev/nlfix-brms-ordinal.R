# What brms 2.23.0 builds for an nl body under each ordinal family, with
# and without thres(), cs() and disc. Prints, per spelling, whether
# make_stancode() succeeds and the lines that compute mu and read the
# thresholds. No compile here; the L1 runs (dev/nlfix-l1-*.R) compile.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
src <- new.env()
sys.source("dev/equiv-specs.R", envir = src)
d <- src$eq_dat("ordinal", seed = 49)
d$y0 <- ordered(ifelse(runif(nrow(d)) < 0.2, 0L, as.integer(d$y)),
                levels = 0:4)
fams <- c("cumulative", "sratio", "cratio", "acat", "hurdle_cumulative")
forms <- list(
  plain = "bf(Y ~ b1 * x, b1 ~ 1, nl = TRUE)",
  b0 = "bf(Y ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE)",
  thres = "bf(Y | thres(3) ~ b1 * x, b1 ~ 1, nl = TRUE)",
  disc = "bf(Y ~ b1 * x, b1 ~ 1, disc ~ 0 + z, nl = TRUE)",
  cs_nlpar = "bf(Y ~ b1 * x, b1 ~ 1 + cs(z), nl = TRUE)",
  cs_body = "bf(Y ~ b1 * x + cs(z), b1 ~ 1, nl = TRUE)",
  nlpar_x = "bf(Y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ 1 + z, nl = TRUE)",
  re = "bf(Y ~ b1 * x, b1 ~ 1 + (1 | g), nl = TRUE)",
  nlf = "bf(Y ~ eta) + nlf(eta ~ b1 * x) + lf(b1 ~ 1)"
)
for (fam in fams) {
  for (nm in names(forms)) {
    yv <- if (fam == "hurdle_cumulative") "y0" else "y"
    txt <- gsub("Y", yv, forms[[nm]], fixed = TRUE)
    if (fam == "hurdle_cumulative" && nm == "thres") next
    if (fam == "hurdle_cumulative") {
      txt <- sub(")$", ", hu ~ 1)", txt)
    }
    f <- eval(parse(text = txt))
    r <- tryCatch(suppressWarnings(suppressMessages(
      stancode(f, data = d, family = get(fam)()))),
      error = function(e) e)
    cat(sprintf("\n== %s / %s: %s\n", fam, nm, txt))
    if (inherits(r, "error")) {
      cat("  REFUSES:", gsub("\n", " ", conditionMessage(r)), "\n")
      next
    }
    ln <- strsplit(r, "\n")[[1]]
    keep <- grep(paste0("nlp_|mu\\[n\\] =|mu =|Intercept|disc|bcs|mucs|",
                        "_lpmf|target \\+="), ln, value = TRUE)
    keep <- unique(trimws(keep))
    cat(paste0("  ", head(keep, 30)), sep = "\n")
  }
}
