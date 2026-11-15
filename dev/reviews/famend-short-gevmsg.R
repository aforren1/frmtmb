# Short check of lane famend, F4: the GEV xi message's confint() call
# runs as written, for a scalar xi, for xi ~ x, and for xi held constant.
#   Rscript dev/reviews/famend-short-gevmsg.R <lib>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
rg <- function(n, mu, s, xi) mu + s * ((-log(stats::runif(n)))^(-xi) - 1) / xi
run <- function(label, fo, d) {
  msgs <- list()
  warns <- character()
  f <- tryCatch(withCallingHandlers(
    frm(fo, family = gen_extreme_value(), data = d),
    message = function(m) {
      msgs[[length(msgs) + 1L]] <<- m
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(f, "error")) {
    cat(label, ": fit error:", conditionMessage(f), "\n")
    return(invisible())
  }
  hit <- Filter(function(m) inherits(m, "frmtmb_gev_xi_nonregular"), msgs)
  cat(sprintf("%s: code %d, xi range %s, warnings %d, nonregular message %d\n",
              label, f$opt$convergence,
              paste(format(range(frm_linpred(f, dpar = "xi",
                                             type = "response")),
                           digits = 4), collapse = " to "),
              length(warns), length(hit)))
  if (length(warns)) cat("   W:", substr(warns, 1, 150), sep = "\n   W: ")
  if (!length(hit)) return(invisible())
  msg <- conditionMessage(hit[[1L]])
  call <- sub(".*(confint\\(fit, parm = \"[^\"]+\", method = \"profile\"\\)).*",
              "\\1", msg)
  cat("   call in message:", call, "\n")
  fit <- f
  ci <- tryCatch(suppressWarnings(eval(parse(text = call))),
                 error = function(e) e)
  if (inherits(ci, "error")) {
    cat("   call ERROR:", conditionMessage(ci), "\n")
  } else {
    print(ci)
  }
}
set.seed(4)
d1 <- data.frame(y = rg(300, 0, 1, -0.63))
run("scalar xi -0.63", y ~ 1, d1)
set.seed(5)
d2 <- data.frame(x = stats::runif(300, -1, 1))
d2$y <- vapply(seq_len(300), function(i) rg(1, 0, 1, -0.65 + 0.1 * d2$x[i]),
               0)
run("xi ~ x", bf(y ~ 1, xi ~ x), d2)
run("xi held at -0.7", bf(y ~ 1, xi = -0.7), d1)
run("mu ~ x, xi held at -0.7", bf(y ~ x, xi = -0.7), d2)
