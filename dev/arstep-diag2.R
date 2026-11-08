# Lane arstep: where the default descent begins on the review's row, against
# start = list(theta = ) at the same theta. Data as dev/arstep-diag.R.
# Args: arm [case n seed]
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
arm <- args[1]
cs <- if (length(args) > 1) args[2] else "gar"
n <- if (length(args) > 2) as.integer(args[3]) else 160L
seed <- if (length(args) > 3) as.integer(args[4]) else 1L
rev_libs(arm)
ns <- asNamespace("frmtmb")
f <- function(x) 2 * (x > 5) - 1
set.seed(74000 + seed * 31 + n)
x <- sort(round(runif(n, 0, 10), 4))
eta <- f(x)
m <- switch(cs,
  student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)), fml = y ~ gp(x), fam = student()),
  gar = list(d = data.frame(x = x, t = seq_len(n), y = eta + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03))), fml = y ~ gp(x) + ar(time = t), fam = gaussian()))
starts <- list()
trace("optimize_obj", tracer = quote(assign("starts", c(get("starts", globalenv()), list(start_par)), envir = globalenv())), print = FALSE, where = ns)
fitit <- function(start = NULL) { a <- list(m$fml, m$d, family = m$fam); if (!is.null(start)) a$start <- start; do.call(rev_fit, a) }
lab <- function(lbl, r) {
  e <- r$fit$estimates
  cat(sprintf("D2 %s ll %.4f conv %d | start %s | end theta %s betad %s beta %s ac %s\n", lbl, as.numeric(logLik(r$fit)), r$fit$opt$convergence,
              paste(names(starts[[1]]), signif(starts[[1]], 5), sep = "=", collapse = " "),
              paste(signif(e$theta, 5), collapse = ","), signif(e$betad, 5), signif(e$beta, 5), signif(e$thetaac %||% NA, 5)))
}
starts <- list(); r0 <- fitit(); lab("default", r0)
w <- ns$gp_start_working(r0$fit$frame)
starts <- list(); r1 <- fitit(list(theta = w$theta)); lab("theta=working", r1)
st <- starts
p0 <- st[[1]]
cat("D2 template par:", paste(names(p0), signif(p0, 5), collapse = " "), "\n")
cat("D2 par_template:\n"); str(r0$fit$frame$par_template[c("beta", "betad", "theta", "thetaac")])
cat("DONE\n")
