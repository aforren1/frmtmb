# Which reader spellings does this build take? Only so that the order script reads
# them and not its own typos.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
set.seed(1)
d <- data.frame(x = runif(40, 0, 6), z = rnorm(40), g = factor(sample(letters[1:5], 40, TRUE)))
d$y <- sin(d$x) + rnorm(40, sd = 0.2)
fit <- suppressWarnings(frm(y ~ z + gp(x) + (1 | g), d))
tr <- function(lab, e) {
  r <- tryCatch({ v <- e; paste(class(v)[1], paste(dim(v), collapse = "x"), length(v)) },
                error = function(err) paste("ERR", substr(conditionMessage(err), 1, 200)))
  cat(sprintf("%-14s %s\n", lab, r))
}
print(names(fit$opt$par))
tr("prof", confint(fit, method = "profile", parm = "z"))
tr("prof2", confint(fit, method = "profile", parm = 2L))
tr("gmeb", getME(fit, "b"))
tr("allfit", frm_allfit(fit))
tr("hyp", hypothesis(fit, "z > 0"))
tr("upd start", update(fit, start = list(theta = fit$estimates$theta)))
tr("boot", frm_bootstrap(fit, nsim = 2))
print(args(frm_bootstrap))
print(args(frm_allfit))
tr("sand", vcov(fit, type = "HC0"))
print(getAnywhere("vcov.frmtmb_fit")$objs[[1]] |> args())
tr("sandwich", sandwich::sandwich(fit))
ex <- ls("package:frmtmb")
print(ex[grepl("boot|allfit|prof|hypoth|sandw|cluster|influ|simulate|linpred|curve", ex)])
b <- frm_bootstrap(fit, nsim = 2); str(b, max.level = 1)
