# Reviewer of lane optrap, overlap with seflags: lane seflags' degree-7
# raw polynomial with (1 | g) (dev/seflags-deg7diag.R's data, seed 7):
# the optimizer's end on one build, what it says, and the units it used.
#   Rscript dev/optrap-rev-deg7.R <lib or none> [deg] [seed]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
deg <- if (length(a) > 1) as.integer(a[2]) else 7L
sd <- if (length(a) > 2) as.integer(a[3]) else deg
cat("lib", dirname(find.package("frmtmb")), "deg", deg, "seed", sd, "\n")
set.seed(sd)
d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + "),
                       "+ (1 | g)"))
w <- character()
f <- withCallingHandlers(suppressMessages(frm(fo, data = d)),
  warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
m <- suppressWarnings(lmer(fo, data = d, REML = FALSE))
cat("code", f$opt$convergence, f$opt$message, "| logLik",
    format(as.numeric(logLik(f)), digits = 10), "lmer",
    format(as.numeric(logLik(m)), digits = 10), "| evals",
    f$opt$evaluations, "\n")
cat("par_units", signif(f$par_units, 3), "\n")
for (x in w) cat("W:", substr(x, 1, 160), "\n")
# the lm() fit of the same fixed part from lmer's fixed effects and the
# frm fit's: which one the frm objective prefers, taken at lmer's point
p <- f$opt$par
nm <- frmtmb:::outer_par_names(f)
p2 <- p
p2[seq_len(deg + 1)] <- fixef(m)
p2[match("sigma_(Intercept)", nm)] <- log(sigma(m))
p2[match("theta_1", nm)] <- log(as.data.frame(VarCorr(m))$sdcor[1])
cat("frm objective at lmer's estimates", format(f$obj$fn(p2), digits = 10),
    "| at its own", format(f$opt$objective, digits = 10), "\n")
