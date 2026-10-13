# Reviewer, re-check (m7): how brms 2.23.0 and the punched lane read a
# `||` term, with and without mo(), with numeric and factor variables.
# brms: the group-level blocks of standata (columns per block, and
# whether the block has correlations); frmtmb: re_blocks, covstruct
# and column names, and logLik where the two spellings should agree.
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(4)
n <- 400
dd <- data.frame(x = sample(0:3, n, TRUE), z = rnorm(n),
                 f = factor(sample(c("a", "b", "c"), n, TRUE)),
                 f2 = factor(sample(c("a", "b"), n, TRUE)),
                 g = factor(sample(paste0("g", 1:25), n, TRUE)))
dd$y <- 1 + 0.5 * dd$x + 0.3 * dd$z + rnorm(25)[dd$g] +
  rnorm(25, 0, 0.5)[dd$g] * dd$z + rnorm(n)
brms_blocks <- function(fo) {
  sc <- brms::stancode(fo, data = dd)
  sd <- brms::standata(fo, data = dd)
  ids <- sort(unique(as.integer(sub("^Z_([0-9]+)_.*", "\\1",
                                    grep("^Z_[0-9]+_", names(sd), value = TRUE)))))
  lines <- strsplit(sc, "\n")[[1]]
  out <- vapply(ids, function(i) {
    m <- sd[[paste0("M_", i)]]
    cor <- any(grepl(paste0("L_", i, "\\b"), lines))
    paste0("block ", i, ": M ", m, if (cor) " correlated" else " uncorrelated")
  }, "")
  pr <- brms::default_prior(fo, data = dd)
  sds <- pr$coef[pr$class == "sd" & nzchar(pr$coef)]
  paste(paste(out, collapse = "; "), "| sd coefs:", paste(sds, collapse = " "))
}
frm_blocks <- function(fo) {
  r <- tryCatch(suppressWarnings(frm(bf(fo), family = gaussian(), data = dd)),
                error = function(e) e)
  if (inherits(r, "error")) return(list(txt = paste("ERROR", conditionMessage(r))))
  txt <- paste(vapply(r$frame$re_blocks, function(b) {
    paste0(b$covstruct, "[", paste(unlist(lapply(b$components, `[[`, "cnms")),
                                   collapse = ","), "]")
  }, ""), collapse = " ; ")
  list(txt = txt, ll = as.numeric(logLik(r)))
}
fos <- list(y ~ z + (1 + z || g), y ~ z + f + (1 + z + f || g),
            y ~ f2 + (f2 || g), y ~ f + (f || g),
            y ~ mo(x) + (1 + mo(x) || g), y ~ mo(x) * f2 + (mo(x) * f2 || g),
            y ~ mo(x) + z + (1 + mo(x) + z || g),
            y ~ mo(x) + f + (1 + mo(x) + f || g))
for (fo in fos) {
  cat("\n==", deparse1(fo), "\n")
  fr <- frm_blocks(fo)
  cat("  frmtmb:", fr$txt, if (!is.null(fr$ll)) paste("| logLik", format(fr$ll, digits = 10)), "\n")
  b <- tryCatch(brms_blocks(fo), error = function(e) paste("ERROR", conditionMessage(e)))
  cat("  brms  :", b, "\n")
}
# the likelihood brms's reading gives for the factor term, written so
# that frmtmb fits it: one diag block over the brms columns
cat("\n== brms's columns written out for frmtmb:\n")
dd$fb <- as.numeric(dd$f == "b"); dd$fc <- as.numeric(dd$f == "c")
for (fo in list(y ~ z + f + (1 + z + fb + fc || g), y ~ z + f + (1 + z + f || g))) {
  fr <- frm_blocks(fo)
  cat("  ", deparse1(fo), ":", fr$txt, "| logLik", format(fr$ll, digits = 10), "\n")
}
