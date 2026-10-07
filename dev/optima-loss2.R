# Lane optima, item 4a: vigport defect 7, brms_multilevel's fit_loss2.
# Is there a better optimum than the default fit reports? The default
# fit, then a long multi-start: 40 starts jittered about the default
# fit's estimates (beta by 10 percent, the covariance parameters by
# N(0, 1)), each refit with nlminb, and the best log-likelihood any of
# them reaches, with its convergence code and max |gradient|.
#   Rscript dev/optima-loss2.R base|lane [nstart]
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
nstart <- if (length(args) > 1) as.integer(args[2]) else 40L
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r6",
  lane = c("C:/Users/adf44/source/r/wt-optima-lib",
           "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
csv <- "dev/optima-log/ClarkTriangle.csv"
if (!file.exists(csv)) {
  utils::download.file(paste0("https://raw.githubusercontent.com/mages/",
                              "diesunddas/master/Data/ClarkTriangle.csv"),
                       csv, quiet = TRUE)
}
loss <- utils::read.csv(csv)
nlform2 <- bf(cum ~ ult * (1 - exp(-(dev / theta)^omega)),
              ult ~ 1 + (1 | ID1 | AY), omega ~ 1 + (1 | ID1 | AY),
              theta ~ 1 + (1 | ID1 | AY), nl = TRUE)
one <- function(start, ctl = frmtmb_control()) {
  w <- character()
  f <- tryCatch(withCallingHandlers(
    frm(nlform2, data = loss, family = gaussian(), start = start,
        control = ctl),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }), error = function(e) conditionMessage(e))
  if (is.character(f)) {
    return(list(ll = NA, code = NA, grad = NA, msg = substr(f, 1, 80),
                fit = NULL))
  }
  list(ll = as.numeric(logLik(f)), code = f$opt$convergence,
       grad = max(abs(f$obj$gr(f$opt$par))),
       msg = paste(f$opt$message, if (length(w)) "| warned"), fit = f,
       evals = f$opt$evals)
}
t0 <- proc.time()[[3]]
d0 <- one(list(beta = c(5000, 1, 45)))
cat(sprintf("default: logLik %.6f code %s max|grad| %.3g (%s) evals %s, %.1f s\n",
            d0$ll, d0$code, d0$grad, d0$msg, d0$evals,
            proc.time()[[3]] - t0))
est <- d0$fit$estimates
set.seed(20261006)
rows <- list()
for (i in seq_len(nstart)) {
  st <- list(beta = est$beta * exp(stats::rnorm(length(est$beta), 0, 0.1)),
             theta = est$theta + stats::rnorm(length(est$theta), 0, 1),
             betad = est$betad + stats::rnorm(length(est$betad), 0, 0.3))
  r <- one(st)
  rows[[i]] <- data.frame(start = i, ll = r$ll, code = r$code,
                          grad = r$grad, msg = r$msg)
}
X <- do.call(rbind, rows)
utils::write.table(X, sprintf("dev/optima-log/loss2-%s.tsv", arm),
                   sep = "\t", quote = FALSE, row.names = FALSE)
cat("starts", nrow(X), "; errors", sum(is.na(X$ll)), "\n")
cat("logLik over starts: max", format(max(X$ll, na.rm = TRUE), digits = 10),
    "; code 0 max", format(max(X$ll[X$code %in% 0], na.rm = TRUE),
                           digits = 10), "\n")
print(table(code = X$code, useNA = "ifany"))
print(utils::head(X[order(-X$ll), ], 10))
cat("default below the best start by",
    format(max(X$ll, na.rm = TRUE) - d0$ll, digits = 4), "\n")
