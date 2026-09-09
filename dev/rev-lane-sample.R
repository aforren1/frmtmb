.libPaths(c("C:/Users/adf44/source/r/rev-lincmt-lib","C:/Users/adf44/source/r/reflib-r2","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb);library(frmtmb.ode);library(frmtmb.sample)})
set.seed(9); ns <- 6L
ts <- c(0.25,0.5,1,2,4,8)
lka <- log(1)+rnorm(ns,0,0.3); lke <- log(0.2)+rnorm(ns,0,0.25)
d <- data.frame(id=factor(rep(seq_len(ns),each=length(ts))), time=rep(ts,ns), dose=100)
i <- as.integer(d$id); ka <- exp(lka[i]); ke <- exp(lke[i])
d$conc <- 100*ka/(10*(ka-ke))*(exp(-ke*d$time)-exp(-ka*d$time)) + rnorm(nrow(d),0,0.2)
f <- bf(conc ~ frm_lincmt(parms=list(ka=exp(lka),ke=exp(lke),V=exp(lV)),
                          times=time, group=id, ncmt=1, depot=TRUE,
                          init=list(depot=dose)),
        lka ~ 1 + (1|id), lke ~ 1, lV ~ 1, nl=TRUE) + gaussian()
t0 <- proc.time()[["elapsed"]]
r <- tryCatch({
  fit <- frm_sample(f, data=d, start=list(beta=c(0, log(0.25), log(8))),
                    chains=1, iter=400, warmup=200, seed=1, refresh=0)
  print(summary(fit))
  "OK"
}, error=function(e) paste("ERROR:", conditionMessage(e)))
cat("\nresult:", r, " seconds", proc.time()[["elapsed"]]-t0, "\n")
