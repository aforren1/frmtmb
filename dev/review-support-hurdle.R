# Reviewer probe (lane support): point-mass families under a
# non-binding trunc(lb = -1)
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
set.seed(404)
n <- 300
y <- rexp(n) * (runif(n) > 0.2)
u <- rbeta(n, 2, 3) * (runif(n) > 0.2)
d <- data.frame(y = y, u = u, lb = -1)
for (cs in list(list("hurdle_gamma", hurdle_gamma(), "y"),
                list("hurdle_lognormal", hurdle_lognormal(), "y"),
                list("zero_inflated_beta", zero_inflated_beta(), "u"))) {
  f0 <- tryCatch(suppressWarnings(frm(as.formula(paste(cs[[3]], "~ 1")), family = cs[[2]], data = d)), error = function(e) conditionMessage(e))
  f1 <- tryCatch(suppressWarnings(frm(as.formula(paste(cs[[3]], "| trunc(lb = lb) ~ 1")), family = cs[[2]], data = d)), error = function(e) conditionMessage(e))
  if (is.character(f0) || is.character(f1)) { cat(cs[[1]], "ERR", substr(if (is.character(f1)) f1 else f0, 1, 150), "\n"); next }
  cat(cs[[1]], "dlogLik", as.numeric(logLik(f1)) - as.numeric(logLik(f0)), " dpar", max(abs(f1$opt$par - f0$opt$par)), "\n")
}
