# Release reviewer of 0.71.0: lane dblbar's || reading inside lane
# famsimplex's simplex and categorical location formulas, against brms
# 2.23.0's standata() Z columns; the location-constant refusal combined
# with ||; and one joint-density check of a fit against brms's Stan
# program (brms_lp_check(), helper-brms.R).
#   Rscript dev/relrev071-dblsx.R [lp]
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = paste0("C:/Users/adf44/source/r/",
                                      "frmtmb-wt-release/dev/relrev071-stan"))
dir.create(Sys.getenv("FRMTMB_STAN_CACHE"), showWarnings = FALSE)
suppressMessages({library(frmtmb); library(testthat)})
cat("lib", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "brms", as.character(packageVersion("brms")), "\n")
tt <- "C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat"
for (h in c("helper-brms.R", "helper-stan.R")) {
  sys.source(file.path(tt, h), envir = globalenv())
}
# brms_block_Z(), frm_block_Z(), same_up_to_levels() from test-double-bar.R
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
eta_b <- 0.3 + 0.5 * dd$x + rnorm(16, 0, 0.5)[dd$g] +
  rnorm(16, 0, 0.3)[dd$g] * dd$x
eta_c <- -0.2 - 0.4 * dd$x
P <- cbind(1, exp(eta_b), exp(eta_c))
P <- P / rowSums(P)
G <- matrix(rgamma(length(P), shape = P * 20), n)
Y <- G / rowSums(G)
colnames(Y) <- c("a", "b", "c")
dd$Y <- Y
dd$yc <- factor(apply(P, 1, function(p) sample(c("a", "b", "c"), 1,
                                               prob = p)))

zcmp <- function(lbl, ff, bf_b, family, family_b) {
  sd_b <- tryCatch(brms_standata(bf_b, data = dd, family = family_b),
                   error = function(e) e)
  fr <- tryCatch(frm(ff, family = family, data = dd, dry_run = "frame"),
                 error = function(e) e)
  if (inherits(sd_b, "error") || inherits(fr, "error")) {
    cat(sprintf("%-46s brms %s | frmtmb %s\n", lbl,
                if (inherits(sd_b, "error")) paste("ERROR",
                  substr(conditionMessage(sd_b), 1, 60)) else "ok",
                if (inherits(fr, "error")) paste("ERROR",
                  substr(conditionMessage(fr), 1, 60)) else "ok"))
    return(invisible())
  }
  ids <- seq_along(grep("^M_", names(sd_b)))
  # blocks of several category predictors share one design, so a block
  # can match several brms blocks: ask for a perfect matching instead
  S <- sapply(ids, function(id) vapply(fr$re_blocks, function(bk) {
    same_up_to_levels(frm_block_Z(fr, bk), brms_block_Z(sd_b, id))
  }, NA))
  S <- matrix(S, length(fr$re_blocks))
  # same component (dpar) on both sides where brms names it
  dp_b <- vapply(ids, function(id) {
    z <- grep(paste0("^Z_", id, "_"), names(sd_b), value = TRUE)[1]
    sub(paste0("^Z_", id, "_(.*)_[0-9]+$"), "\\1", z)
  }, "")
  dp_f <- vapply(fr$re_blocks, function(bk) {
    sub("^.*[.]", "", bk$components[[1L]]$lp_key)
  }, "")
  for (i in seq_len(nrow(S))) for (j in seq_len(ncol(S))) {
    if (grepl("^mu", dp_b[j]) && dp_b[j] != dp_f[i]) S[i, j] <- FALSE
  }
  ok <- nrow(S) == ncol(S) && all(rowSums(S) == 1) && all(colSums(S) == 1)
  nb <- vapply(ids, function(id) sd_b[[paste0("M_", id)]], 1)
  nf <- vapply(fr$re_blocks, function(bk) bk$dim, 1)
  cat(sprintf("%-46s blocks brms %d (dims %s; %s) frmtmb %d (dims %s; %s) | each block = brms's block of its dpar: %s\n",
              lbl, length(ids), paste(nb, collapse = ","),
              paste(dp_b, collapse = ","),
              length(fr$re_blocks), paste(nf, collapse = ","),
              paste(dp_f, collapse = ","), ok))
}
dir_f <- dirichlet()
zcmp("dirichlet mub ~ x + (1 + x || g)",
     bf(Y ~ x, mub ~ x + (1 + x || g)),
     brms::bf(Y ~ x, mub ~ x + (1 + x || g)), dir_f, brms::dirichlet())
zcmp("dirichlet Y ~ x + (1 + f || g) (both mu)",
     bf(Y ~ x + (1 + f || g)), brms::bf(Y ~ x + (1 + f || g)), dir_f,
     brms::dirichlet())
zcmp("dirichlet muc ~ (1 + f || g), refcat b",
     bf(Y ~ x, muc ~ (1 + f || g)), brms::bf(Y ~ x, muc ~ (1 + f || g)),
     dirichlet(refcat = "b"), brms::dirichlet(refcat = "b"))
zcmp("dirichlet mua ~ x + (1 + x | gr(g, cor = FALSE))",
     bf(Y ~ x, mua ~ x + (1 + x | gr(g, cor = FALSE))),
     brms::bf(Y ~ x, mua ~ x + (1 + x | gr(g, cor = FALSE))),
     dirichlet(refcat = "c"), brms::dirichlet(refcat = "c"))
zcmp("logistic_normal Y ~ x + (1 + x || g)",
     bf(Y ~ x + (1 + x || g)), brms::bf(Y ~ x + (1 + x || g)),
     logistic_normal(), brms::logistic_normal())
zcmp("dirichlet2 Y ~ x + (1 + f || g)",
     bf(Y ~ x + (1 + f || g)), brms::bf(Y ~ x + (1 + f || g)),
     brmsfamily("dirichlet2"), brms::brmsfamily("dirichlet2"))
zcmp("categorical mub ~ x + (1 + x || g)",
     bf(yc ~ x, mub ~ x + (1 + x || g)),
     brms::bf(yc ~ x, mub ~ x + (1 + x || g)), categorical(),
     brms::categorical())
zcmp("categorical yc ~ x + (1 + f || g)",
     bf(yc ~ x + (1 + f || g)), brms::bf(yc ~ x + (1 + f || g)),
     categorical(), brms::categorical())

cat("\n-- location constants with ||\n")
ref1 <- function(lbl, expr) {
  r <- tryCatch({ force(expr); "accepted" }, error = function(e)
    paste("REFUSED:", substr(gsub("\n", " ", conditionMessage(e)), 1, 200)))
  cat(sprintf("%-40s %s\n", lbl, r))
}
ref1("dirichlet Y ~ (1 + x || g), muc = 0",
     frm(bf(Y ~ x + (1 + x || g), muc = 0), family = dirichlet(),
         data = dd, dry_run = "frame"))
ref1("dirichlet mub ~ (1 + x || g), muc = 0.5",
     frm(bf(Y ~ x, mub ~ (1 + x || g), muc = 0.5), family = dirichlet(),
         data = dd, dry_run = "frame"))
ref1("categorical mub ~ (1+x||g), muc = 0",
     frm(bf(yc ~ x, mub ~ (1 + x || g), muc = 0), family = categorical(),
         data = dd, dry_run = "frame"))
ref1("dirichlet phi = 30 with || (not a location)",
     frm(bf(Y ~ x + (1 + x || g), phi = 30), family = dirichlet(),
         data = dd, dry_run = "frame"))
cat("brms on bf(Y ~ x + (1 + x || g), muc = 0):",
    tryCatch({
      sc <- brms::stancode(brms::bf(Y ~ x + (1 + x || g), muc = 0),
                           data = dd, family = brms::dirichlet())
      paste("muc constant in Stan code:",
            any(grepl("muc = 0|real muc", strsplit(sc, "\n")[[1]])),
            "| b_muc declared:", any(grepl("vector\\[Kc_muc\\] b_muc|b_muc",
                                          strsplit(sc, "\n")[[1]])))
    }, error = function(e) paste("ERROR", conditionMessage(e))), "\n")

if (identical(commandArgs(TRUE)[1], "lp")) {
  cat("\n-- joint density against brms's Stan program\n")
  ff <- bf(Y ~ x, mub ~ x + (1 + x || g))
  fit <- frm(ff, family = dirichlet(), data = dd,
             control = frmtmb_control(grad_tol = 1e-6, restarts = 3))
  print(VarCorr(fit))
  options(frmtmb.brms_lp_report = TRUE)
  r <- testthat::test_that("lp", {
    brms_lp_check(brms::bf(Y ~ x, mub ~ x + (1 + x || g)), brms::dirichlet(),
                  dd, fit, joint = TRUE)
  })
  cat("lp check passed:", r, "\n")
}
