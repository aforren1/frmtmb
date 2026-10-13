# Reviewer of lane optrap, item 6: the ported suite's fixture 1 (student,
# mo(), s(), arma, sigma ~ Trt) from its default start and from jittered
# starts, with the ledger row brmsfit-methods:326 (fitted sigma all > 0)
# at each end point.
#   Rscript dev/optrap-rev-fix1.R <lib or none> <label> <n starts> [sd]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
lab <- a[2]
nst <- as.integer(a[3])
jsd <- if (length(a) > 3) as.numeric(a[4]) else 0.3
cat(lab, "lib", dirname(find.package("frmtmb")), "\n")
fit0 <- get("brmsfit_example1", envir = asNamespace("brms"))
d <- as.data.frame(fit0$data)
attr(d, "terms") <- NULL
form <- bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) +
             (1 + Trt | visit) + arma(visit, patient), sigma ~ Trt)
one <- function(start = NULL) {
  w <- character()
  f <- tryCatch(withCallingHandlers(
    suppressMessages(frm(form, family = student(), data = d,
                         start = start)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(f, "error")) {
    return(list(row = data.frame(obj = NA, code = NA, min_theta = NA,
                                 h326 = NA, nan_err = NA, nwarn = NA,
                                 warn = substr(conditionMessage(f), 1, 60)),
                fit = NULL))
  }
  fi <- tryCatch(suppressWarnings(fitted(f, dpar = "sigma")),
                 error = function(e) NULL)
  h <- if (is.null(fi)) NA else isTRUE(all(fi > 0))
  list(row = data.frame(obj = f$opt$objective, code = f$opt$convergence,
                        min_theta = min(f$estimates$theta),
                        h326 = h,
                        nan_err = if (is.null(fi)) NA else
                          sum(!is.finite(fi[, 2])),
                        nwarn = length(w),
                        warn = substr(gsub("[\t\n]", " ",
                                           paste(unique(substr(w, 1, 40)),
                                                 collapse = " | ")), 1, 120)),
       fit = f)
}
r0 <- one()
rows <- list(cbind(start = "default", r0$row))
est <- r0$fit$estimates
set.seed(20261008)
for (j in seq_len(nst)) {
  st <- lapply(est, function(v) {
    if (is.numeric(v) && length(v)) v + stats::rnorm(length(v), 0, jsd) else v
  })
  # theta at the default's end point can be -1000 on the log scale; a
  # jittered start there is the same flat edge, so sds start at 0
  if (!is.null(st$theta)) st$theta <- stats::rnorm(length(st$theta), 0, 1)
  rj <- one(st)
  rows[[length(rows) + 1L]] <- cbind(start = paste0("jitter", j), rj$row)
}
X <- do.call(rbind, rows)
X$obj <- sprintf("%.10f", X$obj)
print(X, right = FALSE)
dir.create("dev/optrap-rev-log", showWarnings = FALSE)
write.table(X, sprintf("dev/optrap-rev-log/fix1-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
