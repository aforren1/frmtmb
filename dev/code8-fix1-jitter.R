# Fixture 1 from jittered starts (lane optrap's review, item 6, with
# sd 1): which parameters run away where the objective falls below the
# default stop's 109.34, and the ledger row brmsfit-methods:326 at each.
#   CODE8_ARM=before|after Rscript dev/code8-fix1-jitter.R <label> <n> [sd]
a <- commandArgs(TRUE)
lab <- a[1]; nst <- as.integer(a[2])
jsd <- if (length(a) > 2) as.numeric(a[3]) else 1
source("dev/code8-common.R")
options(width = 200)
fit0 <- get("brmsfit_example1", envir = asNamespace("brms"))
d <- as.data.frame(fit0$data)
attr(d, "terms") <- NULL
form <- bf(count ~ Trt * Age + mo(Exp) + s(Age) + volume + offset(Age) +
             (1 + Trt | visit) + arma(visit, patient), sigma ~ Trt)
one <- function(start = NULL) {
  r <- rfit(form, family = student(), data = d, start = start)
  if (r$err) return(list(row = data.frame(obj = NA, code = NA, h326 = NA,
                                          nan_err = NA, extreme = NA,
                                          conds = substr(conditionMessage(r$f),
                                                         1, 80)), fit = NULL))
  f <- r$f
  fi <- tryCatch(quiet(fitted(f, dpar = "sigma")), error = function(e) NULL)
  p <- f$opt$par
  nm <- fns$outer_par_names(f)
  big <- abs(p) > 8
  list(row = data.frame(obj = f$opt$objective, code = f$opt$convergence,
                        h326 = if (is.null(fi)) NA else isTRUE(all(fi > 0)),
                        nan_err = if (is.null(fi)) NA else
                          sum(!is.finite(fi[, 2])),
                        extreme = paste(sprintf("%s=%.3g", nm[big], p[big]),
                                        collapse = " "),
                        conds = substr(cond_txt(r), 1, 200)),
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
  if (!is.null(st$theta)) st$theta <- stats::rnorm(length(st$theta), 0, 1)
  rows[[length(rows) + 1L]] <- cbind(start = paste0("jitter", j), one(st)$row)
}
X <- do.call(rbind, rows)
write.table(X, sprintf("dev/code8-out/fix1-jitter-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
X$obj <- sprintf("%.6f", X$obj)
print(X[, c("start", "obj", "code", "h326", "nan_err", "extreme")],
      right = FALSE)
cat("h326 TRUE:", sum(X$h326 %in% TRUE), "of", nrow(X), "\n")
cat("DONE\n")
