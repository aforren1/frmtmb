LIB <- Sys.getenv("SKEWINIT_LIB", "C:/Users/adf44/source/r/skewinit-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-skewinit/dev/skewinit-lib.R")
dd <- make_data(1, FALSE)
uf <- frm(bf(y ~ xs, sigma ~ 1, alpha ~ 1), family = skew_normal(),
          data = dd, dry_run = "objective")
print(uf$estimates)
fr <- uf$frame
for (l in fr[["linpreds"]]) {
  cat("resp=", format(l[["resp"]]), " dpar=", l[["dpar"]], " par=",
      l[["par"]], " idx=", paste(l[["idx"]], collapse = ","),
      " Xcols=", paste(colnames(l[["X"]]), collapse = ","),
      " const=", !is.null(l[["constant"]]), "\n")
}
r <- frmtmb:::mu_residuals(fr, fr[["linpreds"]][[1]][["resp"]])
cat("sd(resid)", sd(r), " sd(y)", sd(dd$y), " skew(resid)", skew(r),
    " skew(y)", skew(dd$y), "\n")
cat("names(y):", paste(names(fr[["y"]]), collapse = ","), "\n")
X <- fr[["linpreds"]][[1]][["X"]]
cat("class(X):", paste(class(X), collapse = ","), " is.matrix:", is.matrix(X),
    " dim:", paste(dim(X), collapse = "x"), "\n")
y <- fr[["y"]][["y"]]
cat("class(y):", paste(class(y), collapse = ","), " len", length(y),
    " dim:", paste(dim(y), collapse = "x"), "\n")
z <- try(stats::qr.resid(qr(X), y))
cat("qr.resid class:", paste(class(z), collapse = ","), "\n")
