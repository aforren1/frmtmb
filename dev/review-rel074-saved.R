# Release review of 0.74.0, item 3: fits saved by frmtmb 0.73.0
# (rellib-r11, where cox(), car() and mixture_mvn() were core) read under
# 0.74.0.
#   Rscript dev/review-rel074-saved.R save <dir>            (rellib-r11)
#   Rscript dev/review-rel074-saved.R read <dir> <libs;..> <attach;..>
# <attach> names the extensions to library() first ("none" for none).
a <- commandArgs(trailingOnly = TRUE)
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
if (a[1] == "save") {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r11", user))
  suppressMessages(library(frmtmb))
  cat("save with frmtmb", as.character(packageVersion("frmtmb")), "\n")
  g <- expand.grid(r = 1:4, c = 1:4)
  n <- nrow(g)
  W <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:n) {
    if (abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1) W[i, j] <- 1
  }
  dimnames(W) <- list(paste0("L", 1:n), paste0("L", 1:n))
  set.seed(20261023)
  dd <- data.frame(loc = factor(rep(paste0("L", 1:n), each = 6),
                                levels = paste0("L", 1:n)),
                   x = rnorm(6 * n))
  u <- rnorm(n, 0, 0.7)
  dd$y <- 1 + 0.5 * dd$x + u[as.integer(dd$loc)] + rnorm(nrow(dd), 0, 0.5)
  dd$time <- rexp(nrow(dd), exp(-0.5 + 0.6 * dd$x))
  dd$cen <- rbinom(nrow(dd), 1, 0.2)
  fcar <- frm(bf(y ~ x + car(W, gr = loc, type = "escar")),
              family = gaussian(), data = dd, data2 = list(W = W))
  fcox <- frm(bf(time | cens(cen) ~ x), family = cox(df = 4), data = dd)
  dm <- data.frame(row = seq_len(nrow(datasets::faithful)))
  dm$Y <- data.matrix(datasets::faithful)
  fmvn <- frm(bf(Y ~ 1), family = mixture_mvn(K = 2, D = 2), data = dm)
  for (nm in c("fcar", "fcox", "fmvn")) {
    f <- get(nm)
    cat(nm, "logLik", format(as.numeric(logLik(f)), digits = 12),
        "conv", f$opt$convergence, "\n")
    ref <- list(logLik = as.numeric(logLik(f)), vcov = vcov(f),
                fitted = tryCatch(fitted(f), error = function(e) NULL))
    saveRDS(list(fit = f, data = if (nm == "fmvn") dm else dd, W = W,
                 ref = ref), file.path(a[2], paste0(nm, ".rds")))
  }
  quit(save = "no")
}
.libPaths(c(strsplit(a[3], ";", fixed = TRUE)[[1]], user))
suppressMessages(library(frmtmb))
att <- strsplit(a[4], ";", fixed = TRUE)[[1]]
for (p in setdiff(att, "none")) suppressMessages(library(p, character.only = TRUE))
inst <- function(p) nzchar(system.file(package = p))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "| installed: survival",
    inst("frmtmb.survival"), "spatial", inst("frmtmb.spatial"), "latent",
    inst("frmtmb.latent"), "| attached:", a[4], "\n")
for (nm in c("fcar", "fcox", "fmvn")) {
  s <- readRDS(file.path(a[2], paste0(nm, ".rds")))
  fit <- s$fit
  nd <- s$data[1:3, , drop = FALSE]
  calls <- list(
    print = quote(utils::capture.output(print(fit))),
    summary = quote(utils::capture.output(summary(fit))),
    logLik = quote(stopifnot(isTRUE(all.equal(as.numeric(logLik(fit)),
                                              s$ref$logLik)))),
    vcov = quote(stopifnot(isTRUE(all.equal(vcov(fit), s$ref$vcov)))),
    fitted = quote(if (!is.null(s$ref$fitted))
      stopifnot(isTRUE(all.equal(fitted(fit), s$ref$fitted))) else fitted(fit)),
    predict_new = quote(predict(fit, newdata = nd)),
    simulate = quote(simulate(fit, nsim = 1, seed = 1)),
    VarCorr = quote(VarCorr(fit)),
    update = quote(update(fit, data = s$data)))
  cat("==", nm, "\n")
  for (k in names(calls)) {
    r <- tryCatch({ eval(calls[[k]]); "ok" }, error = function(e) {
      kk <- if (inherits(e, "frmtmb_extension_missing")) "REFUSED[ext]" else
        "ERR"
      paste(kk, substr(gsub("\\s+", " ", conditionMessage(e)), 1, 150))
    }, warning = function(w) paste("WARN", substr(conditionMessage(w), 1, 100)))
    cat(sprintf("  %-12s %s\n", k, r))
  }
}
