# Reviewer, lane dblbar: spellings around the new reading and the new
# refusals, on the base and the lane. Records ok/ERROR, blocks, sd/cor
# names and logLik. Data: dev/dblbar-rev-data.R (seed 20261013).
# Usage: Rscript dev/dblbar-rev-edges.R base|lane
# Output: dev/dblbar-rev-log/edges-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-rev-data.R"))
suppressPackageStartupMessages(library(frmtmb))
sink(file.path(wt, sprintf("dev/dblbar-rev-log/edges-%s.txt", arm)),
     split = TRUE)
cat("lib:", find.package("frmtmb"), "\n")
dd <- dblbar_data()
A <- diag(24) * 0.5 + 0.5
dimnames(A) <- list(levels(dd$g), levels(dd$g))
Q <- solve(A)
cvar <- FALSE
cases <- list(
  sigcor = "bf(y ~ x, sigma ~ (1 + x | gr(g, cor = FALSE)))",
  callg = "bf(y ~ x + (1 + x || factor(gch)))",
  inter = "bf(y ~ x + (1 + x || interaction(g1, g2)))",
  corvar = "bf(y ~ x + (1 + x | gr(g, cor = cvar)))",
  corna = "bf(y ~ x + (1 + x | gr(g, cor = NA)))",
  mo = "bf(y ~ mo(xo) + (1 + mo(xo) || g))",
  two = "bf(y ~ x + f + (1 + f || g) + (0 + x | g))",
  hcross = "bf(y ~ x + (1 + x || g:h))",
  diagcor = "bf(y ~ x + diag(1 + x | gr(g, cor = FALSE)))",
  mmcort = "bf(y ~ x + (1 + x || mm(m1, m2, cor = TRUE)))",
  prec = "bf(y ~ x + (1 + x || gr(g, prec = Q)))",
  precf = "bf(y ~ f + (1 + f || gr(g, prec = Q)))",
  covby = "bf(y ~ x + (1 + x || gr(g, cov = A, by = hg)))",
  nest3n = "bf(y ~ x + (1 + x || g1/g2/g3))",
  mvid = paste0("mvbf(bf(y ~ (1 | p | g)), bf(z ~ (1 + x || g)), ",
                "rescor = FALSE)"),
  one = "bf(y ~ x + (1 || g))",
  zero1 = "bf(y ~ x + (0 + x || g))",
  idcort = "bf(y ~ x + (1 + x | p | gr(g, cor = TRUE)))",
  grcorfstud = "bf(y ~ x + (1 + x | gr(g, cor = FALSE, dist = 'student')))",
  covus = "bf(y ~ x + (1 + x | gr(g, cov = A)))",
  diagcov1 = "bf(y ~ x + diag(1 | gr(g, cov = A)))",
  polycov = "bf(y ~ poly(x, 2) + (1 + poly(x, 2) || gr(g, cov = A)))",
  homdcov = "bf(y ~ x + homdiag(0 + x | gr(g, cov = A)))",
  slash_gr = "bf(y ~ x + (1 + x || gr(g1/g2, cor = FALSE)))",
  lincor = "bf(y ~ x + (1 + x || gr(g, cor = 'no')))"
)
for (nm in names(cases)) {
  fit <- tryCatch(suppressWarnings(suppressMessages(
    frm(eval(str2lang(cases[[nm]])), family = gaussian(), data = dd,
        data2 = list(A = A, Q = Q)))),
    error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-10s %-58s ERROR %s\n", nm, cases[[nm]],
                substr(gsub("\n", " ", conditionMessage(fit)), 1, 120)))
    next
  }
  bl <- vapply(fit$frame$re_blocks, function(b) {
    paste0(b$covstruct, "/", b$dim)
  }, "")
  cat(sprintf("%-10s %-58s ok %s ll %.8f\n   %s\n", nm, cases[[nm]],
              paste(bl, collapse = ","), as.numeric(logLik(fit)),
              paste(grep("^(sd|cor)_", variables(fit), value = TRUE),
                    collapse = " ")))
}
fit <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
r <- tryCatch(identical(fitted(fit, newdata = dd[1:5, ],
                               re_formula = ~ (1 + f | gr(g, cor = FALSE))),
                        fitted(fit, newdata = dd[1:5, ])),
              error = function(e) conditionMessage(e))
cat("re_formula gr(cor = FALSE) same as full:", r, "\n")
u <- tryCatch(update(fit, . ~ . + (1 | g1)), error = function(e) e)
cat("update adds (1 | g1):", if (inherits(u, "error")) conditionMessage(u)
    else paste(grep("^sd_", variables(u), value = TRUE), collapse = " "), "\n")
cat("formula(fit):", deparse1(formula(fit)$formula %||% formula(fit)), "\n")
out <- utils::capture.output(print(fit))
cat(out[grep("Group-Level|~g|Multilevel|sd[(]", out)], sep = "\n")
sink()
