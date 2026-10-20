# Review of lane prilabel: gr(id = ) spellings and a by-term in a
# shared Student-t set, frmtmb lane against brms 2.23.0 (no compile).
# Data as test-student-share.R st_data() (seed 21).
.libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
})
sink(file.path(wt, "dev/prilabel-rev-log/misc.txt"), split = TRUE)
set.seed(21)
n <- 360
g <- factor(sample(sprintf("g%02d", 1:24), n, TRUE))
d <- data.frame(g, x = rnorm(n), z = rnorm(n),
                hg = factor(ifelse(as.integer(g) <= 12L, "u", "v")))
d$y <- 1 + 0.5 * d$x + rnorm(24, 0, 0.8)[d$g] + rnorm(n, 0, 0.6)
try_ <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", conditionMessage(e))
})
for (f in list(y ~ x + (1 | gr(g, id = 1)) + (0 + x | gr(g, id = 1)),
               y ~ x + (1 | gr(g, id = "a", cor = FALSE)) +
                 (0 + x | gr(g, id = "a")),
               y ~ x + (1 | gr(g, id = "a", by = hg)) +
                 (0 + x | gr(g, id = "a", by = hg)))) {
  cat("\n--", deparse1(f), "\n")
  b <- try_({
    p <- as.data.frame(brms::default_prior(f, data = d))
    p[p$class %in% c("sd", "cor"), c("class", "coef", "group")]
  })
  r <- try_({
    p <- as.data.frame(frmtmb::default_prior(frmtmb::bf(f), data = d))
    p[p$class %in% c("sd", "cor"), c("class", "coef", "group")]
  })
  rownames(b) <- NULL
  if (is.data.frame(r)) rownames(r) <- NULL
  cat(" identical:", identical(b, r), "\n")
  print(b)
  if (!identical(b, r)) print(r)
}
cat("\n-- Student-t: a gr(by = ) term beside a plain term on g\n")
fb <- y ~ x + (1 | gr(g, by = hg, dist = "student")) +
  (0 + z | gr(g, dist = "student"))
sc <- try_(brms::stancode(fb, data = d))
if (!is.character(sc) || length(sc) > 1) {
  print(grep("dfm_", strsplit(sc, "\n")[[1]], value = TRUE))
} else cat(sc, "\n")
fr <- try_(frm(frmtmb::bf(fb), data = d))
cat("frmtmb:", if (is.character(fr)) fr else
  paste("t_share", deparse(lapply(fr$frame$re_blocks, `[[`, "t_share")),
        "logLik", format(as.numeric(logLik(fr)), digits = 10)), "\n")
sink()
