# Lane optrap, punch round 2: the review's wide high-R^2 design
# (dev/optrap-rev2-b1.R wide, n = 101000, 50 columns, (1 | g)), seed 1,
# which stops 9107.66 below lmer at code 0 with no warning on 0.69.0.
# Where the fit stops, against lmer; the gradient there in natural
# units; what a restart of nlminb from the stop does.
#   Rscript dev/optrap-p2-wide1.R <lib> [seed]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
s <- if (length(a) > 1) as.integer(a[2]) else 1L
cat("lib", dirname(find.package("frmtmb")), "seed", s, "\n")
set.seed(s)
n <- 101000
p <- 49
Z <- matrix(rnorm(n * p), n, p)
colnames(Z) <- paste0("z", seq_len(p))
d <- data.frame(Z, g = factor(rep(seq_len(1000), length.out = n)))
d$y <- 1 + 1e4 * d$z1 + rnorm(1000, 0, 0.5)[d$g] + rnorm(n)
fo <- as.formula(paste("y ~", paste(colnames(Z), collapse = " + "),
                       "+ (1 | g)"))
w <- character()
t0 <- proc.time()[[3]]
f <- withCallingHandlers(frm(fo, data = d), warning = function(x) {
  w <<- c(w, conditionMessage(x))
  invokeRestart("muffleWarning")
})
cat("fit secs", proc.time()[[3]] - t0, "| warnings", length(w), "\n")
l4 <- lmer(fo, data = d, REML = FALSE)
o <- f$opt
cat("code", o$convergence, o$message, "| iterations", o$iterations,
    "| evaluations", deparse(o$evaluations), "\n")
cat("logLik - lmer", as.numeric(logLik(f)) - as.numeric(logLik(l4)), "\n")
u <- f$par_units %||% rep(1, length(o$par))
g <- drop(f$obj$gr(o$par))
nm <- names(o$par)
cat("max |grad| raw", max(abs(g)), "at", which.max(abs(g)),
    "| natural", max(abs(g * u)), "\n")
b <- fixef_by_dpar(f)$mu
bl <- lme4::fixef(l4)
cat("max |beta - lmer|", max(abs(b - bl)), "at", names(b)[which.max(abs(b - bl))],
    "| z1", b[["z1"]], bl[["z1"]], "\n")
cat("sigma", exp(fixef_by_dpar(f)$sigma), "lmer", sigma(l4),
    "| group sd", exp(o$par[nm == "theta"]),
    "lmer", attr(VarCorr(l4)$g, "stddev"), "\n")
cat("start (par_template):\n")
pt <- par_template(fo, data = d)
print(head(pt, 4))
r <- nlminb(o$par, f$obj$fn, f$obj$gr, scale = 1 / u,
            control = list(iter.max = 1000, eval.max = 1000))
cat("restart from the stop: code", r$convergence, r$message,
    "| objective change", r$objective - o$objective, "| iterations",
    r$iterations, "\n")
