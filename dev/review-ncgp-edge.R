# Item 3: edges of the frame's order and the start. One position, two, all tied, gr = FALSE on all
# tied rows, an NA in x, update() with new data, start = list(b = ) of a fit, kept levels, a start
# where the factor fails, frm_allfit's values on both orders. Each arm: lane and base. Arg: arm.
arm <- commandArgs(TRUE)[1]
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs(arm)
rr <- function(lab, e) {
  w <- character()
  v <- tryCatch(withCallingHandlers(e, warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") }), error = function(err) structure(list(conditionMessage(err)), class = "err"))
  if (inherits(v, "err")) { cat(sprintf("RE2 %-28s ERR %s\n", lab, substr(v[[1]], 1, 170))); return(invisible(NULL)) }
  if (inherits(v, "frmtmb_fit")) {
    th <- round(v$estimates$theta, 3)
    cat(sprintf("RE2 %-28s ll %.4f conv %d w %d theta %s\n", lab, logLik(v), v$opt$convergence, length(w), paste(th, collapse = ",")))
    return(invisible(v))
  }
  cat(sprintf("RE2 %-28s %s w %d\n", lab, paste(utils::head(signif(suppressWarnings(as.numeric(unlist(v))), 5), 7), collapse = " "), length(w)))
  invisible(v)
}
set.seed(1)
d <- data.frame(x = runif(30, 0, 6), z = rnorm(30)); d$y <- sin(d$x) + rnorm(30, sd = 0.2)
rr("one position", frm(y ~ gp(x), transform(d, x = 2)))
rr("one position gr=FALSE", frm(y ~ gp(x, gr = FALSE), transform(d, x = 2)))
rr("two positions", frm(y ~ gp(x), transform(d, x = rep(c(1, 2), 15))))
rr("three, gr=FALSE", frm(y ~ gp(x, gr = FALSE), transform(d, x = rep(c(1, 2, 3), 10))))
dn <- d; dn$x[3] <- NA
rr("NA in x", frm(y ~ gp(x), dn))
f <- rr("base fit", frm(y ~ gp(x), d))
rr("update newdata", update(f, newdata = d[1:20, ]))
rr("update start b (own)", update(f, start = list(b = f$estimates$b)))
rr("start b wrong length", frm(y ~ gp(x), d, start = list(b = rnorm(5))))
rr("start theta huge sd", frm(y ~ gp(x), d, start = list(theta = c(800, 0))))
rr("start theta rho tiny", frm(y ~ gp(x), d, start = list(theta = c(0, -800))))
rr("start theta rho huge", frm(y ~ gp(x), d, start = list(theta = c(0, 800))))
rr("x constant-ish 1e-12", frm(y ~ gp(x), transform(d, x = 1 + 1e-12 * x)))
rr("x large offset 1e9", frm(y ~ gp(x), transform(d, x = 1e9 + x)))
rr("x*1e-150", frm(y ~ gp(x), transform(d, x = 1e-150 * x)))
rr("x*1e150", frm(y ~ gp(x), transform(d, x = 1e150 * x)))
rr("y*1e-150", frm(y ~ gp(x), transform(d, y = 1e-150 * y)))
rr("y*1e150", frm(y ~ gp(x), transform(d, y = 1e150 * y)))
rr("y constant", frm(y ~ gp(x), transform(d, y = 3)))
a <- tryCatch(frm_allfit(f), error = function(e) e)
if (inherits(a, "error")) cat("allfit ERR", conditionMessage(a), "\n") else { cat("allfit names", names(a), "\n"); print(utils::head(a$table %||% a$results %||% a[[1]], 8)) }
rr("simulate field sd", { s <- simulate(f, nsim = 200, seed = 2); apply(as.matrix(s), 1, sd)[1:4] })
rr("frm_simulate", frm_simulate(y ~ gp(x), data = d[1:5, ], parameters = list(beta = 0, betad = 0, theta = c(0, 0)), nsim = 1, seed = 1))
cat("DONE\n")
