# Reviewer: brms 2.23.0 facts for the four claimed defects and the rho
# interval, with no C++ compile (stancode, standata, stanc type-check).
# Rscript dev/review-sarfcor-brms-nocompile.R
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages({library(brms); library(rstan)})
cat("brms", format(packageVersion("brms")), "rstan",
    format(packageVersion("rstan")), "\n")
lattice <- function(k, queen = FALSE) {
  g <- expand.grid(r = 1:k, c = 1:k)
  n <- nrow(g)
  outer(1:n, 1:n, function(i, j) {
    dr <- abs(g$r[i] - g$r[j]); dc <- abs(g$c[i] - g$c[j])
    as.numeric(if (queen) pmax(dr, dc) == 1 else dr + dc == 1)
  })
}
try_ <- function(lab, expr) {
  out <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("---", lab, ":", if (is.character(out)) out else format(out), "\n")
}
set.seed(5)

## defect 1: fcor with an indefinite M
n <- 6
V <- 0.6^abs(outer(1:n, 1:n, "-"))
Vi <- V; Vi[1, 2] <- Vi[2, 1] <- 1.2
cat("indefinite M eigen min:", min(eigen(Vi)$values), "\n")
d <- data.frame(y = rnorm(n), x = rnorm(n))
try_("standata fcor(indefinite)", {
  s <- standata(y ~ x + fcor(Vi), data = d, data2 = list(Vi = Vi))
  paste("accepted, Mfcor", nrow(s$Mfcor), "x", ncol(s$Mfcor))
})
code <- stancode(y ~ x + fcor(Vi), data = d, data2 = list(Vi = Vi))
cat("fcor Stan lines with Mfcor:\n")
cat(grep("Mfcor", strsplit(code, "\n")[[1]], value = TRUE), sep = "\n")

## defect 2: row-standardized W of an irregular graph
for (k in c(3, 8)) {
  A <- lattice(k)
  W <- A / rowSums(A)
  ev <- eigen(W)$values
  dd <- data.frame(y = rnorm(nrow(W)), x = rnorm(nrow(W)))
  s <- standata(y ~ x + sar(W), data = dd, data2 = list(W = W))
  cat(sprintf("rook %dx%d row-std: eigen complex %s, max|Im| %.3g; standata eigenMsar class %s\n",
              k, k, is.complex(ev), if (is.complex(ev)) max(abs(Im(ev))) else 0,
              class(s$eigenMsar)[1]))
  A <- lattice(k, queen = TRUE)
  W <- A / rowSums(A)
  ev <- eigen(W)$values
  cat(sprintf("queen %dx%d row-std: eigen complex %s, range(Re) %s; brms interval (%.4f, %.4f), 1/ev interval (%.4f, %.4f)\n",
              k, k, is.complex(ev), paste(format(range(Re(ev)), digits = 4), collapse = " "),
              min(Re(ev)), max(Re(ev)), 1 / min(Re(ev)), 1 / max(Re(ev))))
}

## defect 3: unstandardized W
A <- lattice(4)
dd <- data.frame(y = rnorm(16), x = rnorm(16))
s <- standata(y ~ x + sar(A), data = dd, data2 = list(A = A))
cat(sprintf("0/1 adjacency 4x4: brms eigen range (%.4f, %.4f); valid (1/min, 1/max) = (%.4f, %.4f)\n",
            min(s$eigenMsar), max(s$eigenMsar), 1 / min(s$eigenMsar),
            1 / max(s$eigenMsar)))
code <- stancode(y ~ x + sar(A), data = dd, data2 = list(A = A))
cat(grep("lagsar;|log1m|eigenMsar", strsplit(code, "\n")[[1]], value = TRUE),
    sep = "\n")

## defect 4: sar with addition terms: does brms's code type-check?
A <- lattice(4); W <- A / rowSums(A)
dd$w <- runif(16, 0.5, 1.5); dd$se <- 0.3; dd$cen <- 0L; dd$lb <- -10
for (f in list(y | weights(w) ~ x + sar(W), y | se(se, sigma = TRUE) ~ x + sar(W),
               y | cens(cen) ~ x + sar(W), y | trunc(lb = lb) ~ x + sar(W),
               y ~ x + sar(W))) {
  lab <- deparse1(f)
  code <- tryCatch(stancode(f, data = dd, data2 = list(W = W)),
                   error = function(e) paste("STANCODE ERROR", conditionMessage(e)))
  if (startsWith(code, "STANCODE ERROR")) { cat(lab, ":", code, "\n"); next }
  r <- tryCatch({stanc(model_code = code, verbose = FALSE); "type-checks"},
                error = function(e) paste("stanc ERROR:",
                                          substr(gsub("\n", " ", conditionMessage(e)), 1, 220)))
  cat(lab, ":", r, "\n")
}

## prior: default_prior for lagsar / errorsar
print(default_prior(y ~ x + sar(W), data = dd, data2 = list(W = W)))
print(default_prior(y ~ x + sar(W, type = "error"), data = dd,
                    data2 = list(W = W)))

## new data with the SAME number of rows as M: a counterfactual
fs <- brm(y ~ x + sar(W), data = dd, data2 = list(W = W), empty = TRUE)
nd <- dd; nd$x <- nd$x + 1
try_("sar newdata of the same 16 rows, no newdata2",
     paste(dim(standata(fs, newdata = nd)$Msar), collapse = " x "))
ff <- brm(y ~ x + fcor(V16), data = dd,
          data2 = list(V16 = 0.5^abs(outer(1:16, 1:16, "-"))), empty = TRUE)
try_("fcor newdata of the same 16 rows, no newdata2",
     paste(dim(standata(ff, newdata = nd)$Mfcor), collapse = " x "))
