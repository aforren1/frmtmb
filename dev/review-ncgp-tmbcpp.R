# Item 3: the same model as a plain TMB C++ template, no RTMB, no frmtmb. Seed 20261102.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
cat("lib:", dirname(system.file(package = "TMB")), "TMB", as.character(packageVersion("TMB")), "\n")
set.seed(20261102)
n0 <- 60L
x0 <- round(runif(n0, 0, 6), 3); e0 <- rnorm(n0)
y <- (sin(2 * x0) + 1e-2 * e0)[1:40]; x <- x0[1:40]; n <- 40L
u <- sort(unique(x)); m <- length(u); J <- match(x, u); D2 <- outer(u, u, "-")^2
Sys.setenv(TMP = "C:/Users/adf44/AppData/Local/Temp", TEMP = "C:/Users/adf44/AppData/Local/Temp", TMPDIR = "C:/Users/adf44/AppData/Local/Temp")
tdir <- "C:/Users/adf44/source/r/rev-ncgp-src/tmbcpp"
dir.create(tdir, recursive = TRUE, showWarnings = FALSE)
cpp <- file.path(tdir, "revdefect.cpp")
writeLines(c(
  "#include <TMB.hpp>",
  "template<class Type> Type objective_function<Type>::operator() () {",
  "  DATA_VECTOR(y); DATA_IVECTOR(J); DATA_MATRIX(D2);",
  "  PARAMETER(mu); PARAMETER(ls); PARAMETER(lsd); PARAMETER(lrho); PARAMETER_VECTOR(b);",
  "  int m = D2.rows(); matrix<Type> K(m, m);",
  "  for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) K(i, j) = exp(-D2(i, j) / (Type(2) * exp(Type(2) * lrho))) + (i == j ? Type(1e-12) : Type(0));",
  "  Eigen::LLT< Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> > llt(K);",
  "  matrix<Type> L = llt.matrixL();",
  "  vector<Type> v = L * b.matrix();",
  "  Type nll = -dnorm(b, Type(0), Type(1), true).sum();",
  "  for (int i = 0; i < y.size(); i++) nll -= dnorm(y(i), mu + exp(lsd) * v(J(i)), exp(ls), true);",
  "  return nll;",
  "}"), cpp)
owd <- setwd(tdir)
TMB::compile("revdefect.cpp", framework = "TMBad")
dyn.load(TMB::dynlib("revdefect"))
at <- function(o) { p <- o$par; p[names(p) == "lsd"] <- 1.08396; p[names(p) == "lrho"] <- 0.36733; p[names(p) == "ls"] <- -4.57404; p }
ng <- function(o, p, h = 1e-4) vapply(seq_along(p), function(i) { e <- replace(numeric(length(p)), i, h); (o$fn(p + e) - o$fn(p - e)) / (2 * h) }, 0)
for (comp in 0:1) {
  TMB::config(tmbad.sparse_hessian_compress = comp, DLL = "revdefect")
  o <- TMB::MakeADFun(list(y = y, J = J - 1L, D2 = D2), list(mu = 0, ls = log(0.01), lsd = 0.5, lrho = 0.3, b = numeric(m)),
                      random = c("mu", "b"), DLL = "revdefect", silent = TRUE)
  force(o$env$spHess)
  p <- at(o)
  g <- o$gr(p); fd <- ng(o, p)
  cat(sprintf("RTC TMB-C++ compress %d fn %.8f gr %s | fd %s\n", comp, o$fn(p), paste(sprintf("%.6f", g), collapse = " "), paste(sprintf("%.6f", fd), collapse = " ")))
}
TMB::config(tmbad.sparse_hessian_compress = 0L, DLL = "revdefect")
setwd(owd)
cat("DONE\n")
