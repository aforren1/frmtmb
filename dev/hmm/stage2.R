# Stage 2: covariate-dependent transitions and the three init modes.
suppressMessages(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-hmm",
                                   quiet = TRUE))
source("C:/Users/adf44/source/r/frmtmb-wt-hmm/dev/hmm/hmm-common.R")
D <- "C:/Users/adf44/source/r/frmtmb-wt-hmm/dev/hmm/"
b1 <- readRDS(paste0(D, "probeB1.rds"))
d1 <- readRDS(paste0(D, "probeD1.rds"))
cat("probeD1:", paste(names(d1), collapse = ", "), "\n")
flat <- function(f) unlist(fixef(f))

dd <- b1$dat
fb <- frm(bf(y ~ 1),
          family = hmm(K = 2, gaussian(), time = t, group = g,
                       init = "estimated", trans = ~x),
          data = dd)
cat("B1 frm logLik      :", sprintf("%.9f", as.numeric(logLik(fb))), "\n")
cat("B1 probe reference :", sprintf("%.9f", b1$ll), "\n")
cat("B1 diff            :", abs(as.numeric(logLik(fb)) - b1$ll), "\n")
cat("df                 :", attr(logLik(fb), "df"), "\n")
print(round(flat(fb), 6))
cat("hmm_ldel:", fb$estimates[["hmm_ldel"]], "\n")

## numeric forward at the estimates, with time-varying transitions
e <- flat(fb)
mu <- c(e[["mu1.(Intercept)"]], e[["mu2.(Intercept)"]])
sg <- exp(c(e[["sigma1.(Intercept)"]], e[["sigma2.(Intercept)"]]))
B <- matrix(0, 2, 2); Bx <- matrix(0, 2, 2)
B[1, 2] <- e[["tr12.(Intercept)"]]; Bx[1, 2] <- e[["tr12.x"]]
B[2, 2] <- e[["tr22.(Intercept)"]]; Bx[2, 2] <- e[["tr22.x"]]
lpm <- lpmat_gauss(dd$y, mu, sg)
dl <- softmax0(fb$estimates[["hmm_ldel"]])
rows <- hmm_seq_index(dd$g, dd$t)
Gof <- function(r) tpm_tv_num(B, Bx, dd$x[r], 2L)
llnum <- sum(vapply(rows, function(r) fwd_num_tv(lpm, Gof, r, dl),
                    numeric(1)))
cat("B1 numeric forward :", sprintf("%.9f", llnum),
    " diff:", abs(llnum - as.numeric(logLik(fb))), "\n\n")

## depmixS4 with transition covariates
if (requireNamespace("depmixS4", quietly = TRUE)) {
  set.seed(4)
  dm <- depmixS4::depmix(y ~ 1, data = dd, nstates = 2, transition = ~x,
                         ntimes = as.integer(table(dd$g)))
  best <- -Inf; bp <- NULL
  for (i in 1:8) {
    ff <- try(suppressMessages(depmixS4::fit(
      dm, verbose = FALSE,
      emcontrol = depmixS4::em.control(random.start = TRUE, tol = 1e-12,
                                       maxit = 5000))), silent = TRUE)
    if (!inherits(ff, "try-error")) {
      v <- as.numeric(depmixS4::logLik(ff))
      if (v > best) { best <- v; bp <- depmixS4::getpars(ff) }
    }
  }
  cat("depmixS4 logLik:", sprintf("%.9f", best),
      " diff:", abs(best - as.numeric(logLik(fb))), "\n")
  cat("depmixS4 pars:", paste(round(bp, 5), collapse = " "), "\n")
  cat("frm tr betas :",
      paste(round(c(B[1, 2], Bx[1, 2], B[2, 2], Bx[2, 2]), 5),
            collapse = " "), "\n\n")
}

## --- init modes on probeD1's data ------------------------------------
dq <- d1$dat
fst <- frm(bf(y ~ 1),
           family = hmm(K = 2, gaussian(), time = t, group = g,
                        init = "stationary"), data = dq)
fes <- frm(bf(y ~ 1),
           family = hmm(K = 2, gaussian(), time = t, group = g,
                        init = "estimated"), data = dq)
fun <- frm(bf(y ~ 1),
           family = hmm(K = 2, gaussian(), time = t, group = g,
                        init = "uniform"), data = dq)
for (nm in c("stationary", "estimated", "uniform")) {
  f <- get(paste0("f", substr(nm, 1, 2)))
  cat(sprintf("%-11s logLik %.8f  df %d\n", nm, as.numeric(logLik(f)),
              attr(logLik(f), "df")))
}
cat("probe F3 reference: free -1608.76549264 (df 7), ",
    "stationary -1609.41007570 (df 6)\n")

## the taped stationary solve against the numeric stationary forward
e <- flat(fst)
mu <- c(e[["mu1.(Intercept)"]], e[["mu2.(Intercept)"]])
sg <- exp(c(e[["sigma1.(Intercept)"]], e[["sigma2.(Intercept)"]]))
lg <- c(e[["tr12.(Intercept)"]], e[["tr22.(Intercept)"]])
G <- rbind(c(1, exp(lg[1])) / (1 + exp(lg[1])),
           c(1, exp(lg[2])) / (1 + exp(lg[2])))
dl <- stat_dist(G)
lpm <- lpmat_gauss(dq$y, mu, sg)
rows <- hmm_seq_index(dq$g, dq$t)
llnum <- sum(vapply(rows, function(r) fwd_num(lpm[r, , drop = FALSE], G, dl),
                    numeric(1)))
cat("stationary: taped", sprintf("%.9f", as.numeric(logLik(fst))),
    " numeric", sprintf("%.9f", llnum),
    " diff", abs(llnum - as.numeric(logLik(fst))), "\n")

## refusal: stationary + transition covariates
r <- try(frm(bf(y ~ 1),
             family = hmm(K = 2, gaussian(), time = t, group = g,
                          init = "stationary", trans = ~x),
             data = dd), silent = TRUE)
cat("\nrefusal:", conditionMessage(attr(r, "condition")), "\n")
