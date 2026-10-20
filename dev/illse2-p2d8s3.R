# Lane illse2, item 2: degree 8 seed 3 of the 80 sd-at-zero fits, whose
# nlminb code after the sd escape differs between BLAS: the warning.
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
set.seed(3); n <- 300; deg <- 8
d <- data.frame(x = stats::runif(n, 1, 2), g = factor(rep(1:20, length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
xn <- paste0("x", seq_len(deg))
for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "), "+ (1 | g)"))
r <- cap(frm(fo, data = d))
cat("code", r$v$opt$convergence, r$v$opt$message, "| escape",
    format(r$v$opt$sd_escape), "| polish", format(r$v$opt$polish), "\n")
for (w in r$w) cat("WARN:", w, "\n")
