# Does the post-fit surface work through frm_lincmt()?
.libPaths(c("C:/Users/adf44/source/r/lincmt-lib","C:/Users/adf44/source/r/reflib-r2","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb);library(frmtmb.ode)})
set.seed(77); ns <- 8L
ts <- c(0.25,0.5,1,2,4,8,12)
lka <- log(1)+rnorm(ns,0,0.3); lke <- log(0.2)+rnorm(ns,0,0.25)
d <- data.frame(id=factor(rep(seq_len(ns),each=length(ts))), time=rep(ts,ns), dose=100)
i <- as.integer(d$id); ka <- exp(lka[i]); ke <- exp(lke[i])
d$conc <- 100*ka/(10*(ka-ke))*(exp(-ke*d$time)-exp(-ka*d$time)) + rnorm(nrow(d),0,0.2)
fit <- frm(bf(conc ~ frm_lincmt(parms=list(ka=exp(lka),ke=exp(lke),V=exp(lV)),
                                times=time, group=id, ncmt=1, depot=TRUE,
                                init=list(depot=dose)),
              lka ~ 1 + (1|id), lke ~ 1 + (1|id), lV ~ 1, nl=TRUE) + gaussian(),
           data=d, start=list(beta=c(0,log(0.25),log(8))), se=TRUE)
chk <- function(nm, e) {
  r <- tryCatch({ v <- e; paste0("OK  len=", length(v), " range=",
     paste(format(range(as.numeric(v)), digits=4), collapse=" ")) },
    error=function(err) paste("ERROR:", conditionMessage(err)))
  cat(sprintf("%-28s %s\n", nm, r))
}
chk("fitted()", fitted(fit))
chk("residuals()", residuals(fit))
chk("predict(), training data", predict(fit))
nd <- data.frame(id=factor(rep(1:2, each=5), levels=levels(d$id)),
                 time=rep(c(0.5,1,3,6,10),2), dose=100)
chk("predict(newdata)", predict(fit, newdata=nd))
chk("simulate()", unlist(simulate(fit, nsim=2)))
chk("simulate() ode-free", { s2 <- simulate(fit, nsim=2); unlist(s2) })
chk("logLik()", logLik(fit))
cat("\ndiagnose():\n"); print(diagnose(fit))
