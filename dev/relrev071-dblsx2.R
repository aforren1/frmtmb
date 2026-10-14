# Release reviewer: why frm_block_Z() and brms_block_Z() disagree on a
# simplex location's || block (helper artifact or real?).
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(testthat)})
tt <- "C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat"
sys.source(file.path(tt, "helper-brms.R"), envir = globalenv())
src <- readLines(file.path(tt, "test-double-bar.R"))
a <- grep("^brms_block_Z <- ", src)
b <- grep("^test_that\\(\"every \\|\\| case has brms", src)
eval(parse(text = src[a:(b - 1L)]), envir = globalenv())
sys.source("C:/Users/adf44/source/r/frmtmb-wt-release/dev/relrev071-zpatch.R", envir = globalenv())
set.seed(11)
n <- 400
dd <- data.frame(x = rnorm(n), z = rnorm(n),
                 g = factor(rep(sprintf("g%02d", 1:16), length.out = n)),
                 f = factor(sample(c("a", "b", "c"), n, TRUE)))
Y <- matrix(rgamma(3 * n, 2), n); Y <- Y / rowSums(Y)
colnames(Y) <- c("a", "b", "c"); dd$Y <- Y
dd$yc <- factor(sample(c("a", "b", "c"), n, TRUE))
cmp <- function(lbl, ff, bfb, fam, famb) {
  sd_b <- brms_standata(bfb, data = dd, family = famb)
  fr <- frm(ff, family = fam, data = dd, dry_run = "frame")
  ids <- seq_along(grep("^M_", names(sd_b)))
  cat("==", lbl, "\n")
  for (bk in fr$re_blocks) {
    lpk <- vapply(bk$components, function(cc) cc$lp_key, "")
    Zf <- frm_block_Z(fr, bk)
    cat("  frmtmb block: comps", paste(lpk, collapse = ","), "dim", bk$dim,
        "levels", bk$n_levels, "Zf", dim(Zf), "\n")
    for (id in ids) {
      Zb <- brms_block_Z(sd_b, id)
      cat("    brms id", id, "Zb", dim(Zb), "same", same_up_to_levels(Zf, Zb),
          "| max|Zf - Zb| (same order)",
          if (identical(dim(Zf), dim(Zb))) signif(max(abs(Zf - Zb)), 3)
          else NA, "\n")
    }
  }
  # brms's raw Z vectors and J for id 1
  cat("  brms names:", grep("^(Z|J|N|M)_1", names(sd_b), value = TRUE), "\n")
}
cmp("dirichlet mub ~ x + (1 + x || g)", bf(Y ~ x, mub ~ x + (1 + x || g)),
    brms::bf(Y ~ x, mub ~ x + (1 + x || g)), dirichlet(), brms::dirichlet())
cmp("gaussian y ~ x + (1 + x || g)", bf(x ~ z + (1 + z || g)),
    brms::bf(x ~ z + (1 + z || g)), gaussian(), gaussian())
cmp("categorical mub ~ x + (1 + x || g)", bf(yc ~ x, mub ~ x + (1 + x || g)),
    brms::bf(yc ~ x, mub ~ x + (1 + x || g)), categorical(),
    brms::categorical())
cmp("categorical mub ~ x + (1 + x | g)", bf(yc ~ x, mub ~ x + (1 + x | g)),
    brms::bf(yc ~ x, mub ~ x + (1 + x | g)), categorical(),
    brms::categorical())
