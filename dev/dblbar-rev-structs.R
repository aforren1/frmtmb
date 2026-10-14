# Reviewer, lane dblbar: do the correlated structures keep their names,
# summary rows and fits? Records variables(), summary()$random rows,
# VarCorr() element names, confint_varcorr() rows and logLik per model.
# Data: dev/dblbar-rev-data.R (seed 20261013) plus tim, f2, loc.
# Usage: Rscript dev/dblbar-rev-structs.R base|lane
# Output: dev/dblbar-rev-log/structs-<arm>.rds
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-rev-data.R"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
dd <- dblbar_data()
set.seed(99)
dd$tim <- factor(sample(sprintf("t%d", 1:4), nrow(dd), TRUE))
dd$f2 <- factor(sample(c("a", "b"), nrow(dd), TRUE))
dd$loc <- factor(sprintf("L%d", (as.integer(dd$g) - 1L) %% 6L + 1L))
W <- matrix(0, 6, 6)
for (i in 1:5) W[i, i + 1] <- W[i + 1, i] <- 1
W[1, 6] <- W[6, 1] <- 1
dimnames(W) <- list(sprintf("L%d", 1:6), sprintf("L%d", 1:6))
V <- matrix(c(1.2, 0.5, 0.5, 0.8), 2, 2)
A <- diag(24) * 0.5 + 0.5
dimnames(A) <- list(levels(dd$g), levels(dd$g))
models <- list(
  us = bf(y ~ x + (1 + x | g)),
  us3 = bf(y ~ f + (1 + f | g)),
  cs = bf(y ~ 1 + cs(tim + 0 | g)),
  homcs = bf(y ~ 1 + homcs(tim + 0 | g)),
  ar1 = bf(y ~ 1 + ar1(tim + 0 | g)),
  hetar1 = bf(y ~ 1 + hetar1(tim + 0 | g)),
  toep = bf(y ~ 1 + toep(tim + 0 | g)),
  homdiag = bf(y ~ 1 + homdiag(tim + 0 | g)),
  diag = bf(y ~ 1 + diag(tim + 0 | g)),
  rr = bf(y ~ 1 + rr(tim + 0 | g, d = 1)),
  equalto = bf(y ~ 1 + equalto(f2 + 0 | g, V)),
  stud = bf(y ~ x + (1 + x | gr(g, dist = "student"))),
  studcs = bf(y ~ x + cs(x | gr(g, dist = "student"))),
  studdb = bf(y ~ x + (1 + x || gr(g, dist = "student"))),
  studdbf = bf(y ~ f + (1 + f || gr(g, dist = "student"))),
  car = bf(y ~ x + car(W, gr = loc, type = "icar")),
  grby = bf(y ~ x + (1 + x | gr(g, by = hg))),
  mm = bf(y ~ x + (1 + x | mm(m1, m2))),
  id = bf(y ~ x + (1 + x | p | g), sigma ~ (1 | p | g)),
  grcov = bf(y ~ x + (1 + x | gr(g, cov = A))),
  grcovdb = bf(y ~ x + (1 + x || gr(g, cov = A))),
  idcross = bf(y ~ x + f + (1 + x | p | g) + (1 + f | p | g1))
)
res <- list()
for (nm in names(models)) {
  fit <- tryCatch(suppressWarnings(suppressMessages(
    frm(models[[nm]], family = gaussian(), data = dd,
        data2 = list(W = W, V = V, A = A)))),
    error = function(e) e)
  if (inherits(fit, "error")) {
    res[[nm]] <- list(error = conditionMessage(fit))
    cat(nm, "ERROR", conditionMessage(fit), "\n")
    next
  }
  g <- function(expr) tryCatch(expr, error = function(e) {
    paste("ERROR:", conditionMessage(e))
  })
  cv <- g(suppressWarnings(confint_varcorr(fit)))
  res[[nm]] <- list(
    vars = variables(fit),
    blocks = vapply(fit$frame$re_blocks, function(b) {
      paste0(b$covstruct, "/", b$dim)
    }, ""),
    rand = g(lapply(summary(fit)$random, rownames)),
    vc = g(lapply(VarCorr(fit), names)),
    cv = if (is.data.frame(cv)) {
      paste(cv$type, cv$group %||% "", cv$term %||% cv$name %||% "",
            sep = "|")
    } else cv,
    cvest = if (is.data.frame(cv)) cv$estimate %||% cv[[grep(
      "est", names(cv), ignore.case = TRUE)[1]]] else NULL,
    ll = as.numeric(logLik(fit)),
    hyp = g(names(unlist(frmtmb:::hyp_env_vals(
      fit, frmtmb:::hyp_vals_only(fit)$vals,
      frmtmb:::hyp_vals_only(fit)$comp)))))
  cat(sprintf("%-8s %-22s ll %.10f  sd/cor: %s\n", nm,
              paste(res[[nm]]$blocks, collapse = ","), res[[nm]]$ll,
              paste(grep("^(sd|cor)_", res[[nm]]$vars, value = TRUE),
                    collapse = " ")))
}
saveRDS(res, file.path(wt, sprintf("dev/dblbar-rev-log/structs-%s.rds", arm)))
