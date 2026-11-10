# Reviewer probe (lane support): nested quadrature of a GEV at seed 2,
# xi 0.5 (dev/review-support-is.R), with the rounds shown, and the
# offsets prime() sets on doubles at the Laplace optimum.
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
rgev <- function(n, mu, sigma, xi) mu + sigma * ((-log(runif(n)))^(-xi) - 1) / xi
seed <- as.integer(Sys.getenv("NQ_SEED", "2"))
xi <- as.numeric(Sys.getenv("NQ_XI", "0.5"))
set.seed(seed)
ng <- 15
m <- 8
g <- rep(seq_len(ng), each = m)
u <- rnorm(ng, 0, 0.7)
d <- data.frame(g = factor(g), y = rgev(ng * m, 1 + u[g], 1, xi))
d$h <- factor(paste(d$g, rep(1:2, length.out = nrow(d))))
# trace the offsets each prime() leaves, on doubles
ns <- asNamespace("frmtmb")
orig <- get("quad_nest_objective", ns)
trace_obj <- function(...) {
  qo <- orig(...)
  p0 <- qo$prime
  qo$prime <- function(pars) {
    p0(pars)
    st <- environment(p0)$st
    cat(sprintf("  prime: off_i non-finite %d of %d, off_o non-finite %d of %d\n",
                sum(!is.finite(st$off_i)), length(st$off_i),
                sum(!is.finite(st$off_o)), length(st$off_o)))
    invisible(NULL)
  }
  qo
}
unlockBinding("quad_nest_objective", ns)
assign("quad_nest_objective", trace_obj, ns)
r <- tryCatch(suppressWarnings(
  frm(y ~ 1 + (1 | g) + (1 | h), family = gen_extreme_value(), data = d,
      quadrature = TRUE, verbose = TRUE)),
  error = function(e) e)
if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n") else
  print(c(logLik = as.numeric(logLik(r)), r$opt$par), digits = 12)
