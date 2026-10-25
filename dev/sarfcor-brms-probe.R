# What brms 2.23.0 does with sar() and fcor(): families, Stan code,
# data, default priors, the checks on M. No compile.
# Run: Rscript dev/sarfcor-brms-probe.R > dev/sarfcor-out/brms-probe.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
n <- 9
# 3 x 3 rook lattice
xy <- expand.grid(i = 1:3, j = 1:3)
A <- as.matrix(dist(xy)) == 1
A <- A * 1
W <- A / rowSums(A)                 # row-standardized
d <- data.frame(y = rnorm(n), x = rnorm(n), g = factor(rep(1:3, 3)),
                cnt = rpois(n, 3), yb = rbinom(n, 1, 0.5))
V <- 0.5^abs(outer(1:n, 1:n, "-"))  # an AR(1) correlation as fixed M
V2 <- 2 * V                          # a covariance, not a correlation
dd2 <- list(W = W, A = A, V = V, V2 = V2)

try_ <- function(lab, expr) {
  out <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)),
                  warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat("---", lab, "\n")
  if (is.character(out)) cat(out, sep = "\n") else print(out)
  cat("\n")
}

sc <- function(f, fam = gaussian(), data = d, data2 = dd2) {
  stancode(f, data = data, family = fam, data2 = data2)
}
grepl_code <- function(code, pat) {
  l <- strsplit(code, "\n")[[1]]
  l[grepl(pat, l)]
}

cat("==== stancode, gaussian, lag\n")
code <- sc(y ~ x + sar(W, type = "lag"))
print(code)
cat("==== stancode, gaussian, error (likelihood lines)\n")
code <- sc(y ~ x + sar(W, type = "error"))
cat(grepl_code(code, "sar|target|eigen"), sep = "\n")
cat("==== stancode, student, lag (likelihood lines)\n")
code <- sc(y ~ x + sar(W), fam = student())
cat(grepl_code(code, "sar|target|eigen|nu"), sep = "\n")
cat("==== default type\n")
st <- standata(y ~ x + sar(W), data = d, data2 = dd2)
cat("names:", names(st), "\n")
cat("eigenMsar:", format(st$eigenMsar, digits = 6), "\n")
cat("eigen(W) direct:", format(eigen(W)$values, digits = 6), "\n")
cat("Msar equals W:", isTRUE(all.equal(unname(st$Msar), unname(W))), "\n")

cat("==== default_prior lag\n")
print(default_prior(y ~ x + sar(W, type = "lag"), data = d, data2 = dd2))
cat("==== default_prior error\n")
print(default_prior(y ~ x + sar(W, type = "error"), data = d, data2 = dd2))

cat("==== families for sar\n")
fams <- list(gaussian = gaussian(), student = student(), poisson = poisson(),
             lognormal = lognormal(), skew_normal = skew_normal(),
             bernoulli = bernoulli(), Gamma = Gamma(link = "log"),
             exgaussian = exgaussian(), asym_laplace = asym_laplace())
resp_for <- function(nm) switch(nm, poisson = "cnt", bernoulli = "yb",
                                lognormal = , Gamma = "cnt1", "y")
d$cnt1 <- d$cnt + 1
for (nm in names(fams)) {
  f <- as.formula(paste(resp_for(nm), "~ x + sar(W)"))
  try_(paste("sar family", nm), {
    code <- sc(f, fam = fams[[nm]]); "stancode OK"
  })
}
for (nm in names(fams)) {
  f <- as.formula(paste(resp_for(nm), "~ x + fcor(V)"))
  try_(paste("fcor family", nm), {
    code <- sc(f, fam = fams[[nm]]); "stancode OK"
  })
}

cat("==== the checks on M for sar\n")
try_("sar M binary, not row-standardized", {
  st <- standata(y ~ x + sar(A), data = d, data2 = dd2)
  list(eig = range(st$eigenMsar), M_rowsums = rowSums(st$Msar))
})
try_("sar M binary: default_prior bounds",
     default_prior(y ~ x + sar(A), data = d, data2 = dd2))
Wns <- W; Wns[1, 2] <- Wns[1, 2] + 0.3  # not symmetric, not row-standardized
try_("sar M asymmetric", {
  st <- standata(y ~ x + sar(Wns), data = d, data2 = list(Wns = Wns))
  range(st$eigenMsar)
})
# a matrix with complex eigenvalues: a directed 3-cycle on three nodes
Cc <- matrix(0, n, n)
Cc[cbind(1:n, c(2:n, 1))] <- 1
try_("sar M complex eigenvalues (directed cycle)", {
  st <- standata(y ~ x + sar(Cc), data = d, data2 = list(Cc = Cc))
  head(st$eigenMsar)
})
try_("sar M wrong dimension", standata(y ~ x + sar(Wb), data = d,
                                       data2 = list(Wb = diag(3))))
try_("sar M a data frame", standata(y ~ x + sar(Wd), data = d,
                                    data2 = list(Wd = as.data.frame(W))))
try_("sar M a vector", standata(y ~ x + sar(Wv), data = d,
                                data2 = list(Wv = 1:9)))
try_("sar M nonzero diagonal", {
  Wdg <- W + diag(0.1, n)
  st <- standata(y ~ x + sar(Wdg), data = d, data2 = list(Wdg = Wdg))
  range(st$eigenMsar)
})
try_("sar type invalid", standata(y ~ x + sar(W, type = "mixed"), data = d,
                                  data2 = dd2))
try_("sar sigma predicted", sc(bf(y ~ x + sar(W), sigma ~ x)))
try_("sar nu predicted", sc(bf(y ~ x + sar(W), nu ~ x), fam = student()))
try_("sar two terms", sc(y ~ x + sar(W) + sar(W, type = "error")))
try_("sar + fcor", sc(y ~ x + sar(W) + fcor(V)))
try_("sar + ar(cov = FALSE)", {
  d$t <- 1:9; code <- sc(y ~ x + sar(W) + ar(t)); grepl_code(code, "target")
})
try_("sar + ar(cov = TRUE)", {
  code <- sc(y ~ x + sar(W) + ar(t, cov = TRUE)); grepl_code(code, "target")
})
try_("sar + random effect", {
  code <- sc(y ~ x + sar(W) + (1 | g)); grepl_code(code, "target|r_1")
})
try_("sar + weights", sc(y | weights(cnt1) ~ x + sar(W)))
try_("sar + se", sc(y | se(0.1) ~ x + sar(W)))
try_("sar + cens", sc(y | cens(yb) ~ x + sar(W)))
try_("sar + trunc", sc(y | trunc(lb = -10) ~ x + sar(W)))
try_("sar + mi", sc(y | mi() ~ x + sar(W)))
try_("sar + rescor", sc(bf(mvbind(y, x) ~ sar(W)) + set_rescor(TRUE)))
try_("sar in nl formula", sc(bf(y ~ a * x, a ~ 1 + sar(W), nl = TRUE)))
try_("sar in sigma formula", sc(bf(y ~ x, sigma ~ sar(W))))
try_("sar mixture", sc(y ~ x + sar(W), fam = mixture(gaussian, gaussian)))
try_("sar NA rows", {
  dn <- d; dn$x[2] <- NA
  st <- standata(y ~ x + sar(W), data = dn, data2 = dd2)
  dim(st$Msar)
})
try_("sar M with dimnames vs data order", {
  Wn <- W; dimnames(Wn) <- list(letters[1:9], letters[1:9])
  st <- standata(y ~ x + sar(Wn), data = d, data2 = list(Wn = Wn))
  dimnames(st$Msar)
})
try_("sar M sparse Matrix", {
  st <- standata(y ~ x + sar(Ws), data = d,
                 data2 = list(Ws = Matrix::Matrix(W, sparse = TRUE)))
  class(st$Msar)
})

cat("==== fcor\n")
code <- sc(y ~ x + fcor(V))
print(code)
try_("fcor default_prior", default_prior(y ~ x + fcor(V), data = d,
                                         data2 = dd2))
try_("fcor sigma predicted (het)", {
  code <- sc(bf(y ~ x + fcor(V), sigma ~ x)); grepl_code(code, "target|fcor")
})
try_("fcor student", {
  code <- sc(y ~ x + fcor(V), fam = student()); grepl_code(code, "target|fcor")
})
try_("fcor student nu predicted", {
  code <- sc(bf(y ~ x + fcor(V), nu ~ x), fam = student())
  grepl_code(code, "target|fcor")
})
try_("fcor covariance M (2 V)", {
  st <- standata(y ~ x + fcor(V2), data = d, data2 = dd2)
  diag(st$Mfcor)
})
try_("fcor vector M", {
  st <- standata(y ~ x + fcor(v), data = d, data2 = list(v = 1:9 / 9))
  st$Mfcor
})
try_("fcor asymmetric M", {
  Va <- V; Va[1, 2] <- 0.9
  standata(y ~ x + fcor(Va), data = d, data2 = list(Va = Va))
})
try_("fcor indefinite M (one negative eigenvalue)", {
  Vi <- V; Vi[1, 2] <- Vi[2, 1] <- 1.2
  ev <- eigen(Vi)$values
  st <- standata(y ~ x + fcor(Vi), data = d, data2 = list(Vi = Vi))
  list(eigen = ev, accepted = !is.null(st$Mfcor))
})
try_("fcor negative definite M", {
  Vn <- -V
  standata(y ~ x + fcor(Vn), data = d, data2 = list(Vn = Vn))
})
try_("fcor wrong dimension", standata(y ~ x + fcor(Vb), data = d,
                                      data2 = list(Vb = diag(3))))
try_("fcor NA rows", {
  dn <- d; dn$x[2] <- NA
  st <- standata(y ~ x + fcor(V), data = dn, data2 = dd2)
  dim(st$Mfcor)
})
try_("fcor + ar(cov = FALSE)", {
  code <- sc(y ~ x + fcor(V) + ar(t)); grepl_code(code, "target")
})
try_("fcor + cosy", sc(y ~ x + fcor(V) + cosy(t)))
try_("fcor + unstr", sc(y ~ x + fcor(V) + unstr(t, g)))
try_("fcor + weights", sc(y | weights(cnt1) ~ x + fcor(V)))
try_("fcor + se", {
  code <- sc(y | se(0.1) ~ x + fcor(V)); grepl_code(code, "target|se")
})
try_("fcor + mi", sc(y | mi() ~ x + fcor(V)))
try_("fcor + cens", sc(y | cens(yb) ~ x + fcor(V)))
try_("fcor + random effect", {
  code <- sc(y ~ x + fcor(V) + (1 | g)); grepl_code(code, "target")
})
try_("fcor in sigma formula", sc(bf(y ~ x, sigma ~ fcor(V))))
try_("fcor in nl formula", sc(bf(y ~ a * x, a ~ 1 + fcor(V), nl = TRUE)))
try_("fcor mixture", sc(y ~ x + fcor(V), fam = mixture(gaussian, gaussian)))
try_("fcor rescor", sc(bf(mvbind(y, x) ~ fcor(V)) + set_rescor(TRUE)))
try_("fcor mv without rescor", {
  code <- sc(bf(mvbind(y, x) ~ fcor(V)) + set_rescor(FALSE))
  grepl_code(code, "target")
})
try_("fcor M named but not matched to rows", {
  Vn <- V; dimnames(Vn) <- list(letters[9:1], letters[9:1])
  st <- standata(y ~ x + fcor(Vn), data = d, data2 = list(Vn = Vn))
  dimnames(st$Mfcor)
})

cat("==== an unfitted brmsfit and new data\n")
fit_sar <- brm(y ~ x + sar(W), data = d, data2 = dd2, empty = TRUE)
fit_fc <- brm(y ~ x + fcor(V), data = d, data2 = dd2, empty = TRUE)
nd <- d[1:4, ]
try_("sar standata(newdata) with 4 rows",
     dim(standata(fit_sar, newdata = nd)$Msar))
try_("sar standata(newdata, data2 = W[1:4, 1:4])",
     dim(standata(fit_sar, newdata = nd,
                  data2 = list(W = W[1:4, 1:4]))$Msar))
try_("fcor standata(newdata) with 4 rows",
     dim(standata(fit_fc, newdata = nd)$Mfcor))
try_("fcor standata(newdata, data2 = V[1:4, 1:4])",
     dim(standata(fit_fc, newdata = nd,
                  data2 = list(V = V[1:4, 1:4]))$Mfcor))
try_("sar standata(newdata) same 9 rows",
     dim(standata(fit_sar, newdata = d)$Msar))
cat("variables of the empty fits: none sampled; parameter names from stancode\n")
