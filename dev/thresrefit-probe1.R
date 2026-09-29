## Reachability probe on the REFERENCE build (0.64.0).
## Defect A: does a simulate-and-refit path recount ordinal thresholds?
## Defect B: does draw_prior_pars() recycle one draw into every
## threshold of an unordered threshold vector?
##
## Run: Rscript dev/thresrefit-probe1.R
lib <- Sys.getenv("FRMTMB_PROBE_LIB",
                  "C:/Users/adf44/source/r/rellib-r3")
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("frmtmb", as.character(packageVersion("frmtmb")), "from", lib, "\n")

say <- function(...) cat(..., "\n", sep = "")

## ------------------------------------------------------------------
## A data set whose top ordinal category is rare: n = 60, 4 categories,
## the 4th reached by a couple of rows only.
## ------------------------------------------------------------------
mk_dat <- function(seed, n = 60, tau = c(-0.8, 0.4, 2.6), slope = 0.8) {
  set.seed(seed)
  x <- rnorm(n)
  eta <- slope * x
  p <- cbind(plogis(tau[1] - eta),
             plogis(tau[2] - eta) - plogis(tau[1] - eta),
             plogis(tau[3] - eta) - plogis(tau[2] - eta),
             1 - plogis(tau[3] - eta))
  y <- apply(p, 1, function(pr) sample.int(4L, 1L, prob = pr))
  data.frame(x = x, y = y)
}

dd <- mk_dat(101)
say("A: data seed 101, table(y) = ", paste(table(dd$y), collapse = "/"))

fit <- frm(bf(y ~ x), family = cumulative(), data = dd)
say("A: fitted thresholds = ", length(fit$estimates$tau_raw))
say("A: fixef names = ", paste(names(fixef(fit, flatten = TRUE)),
                               collapse = ", "))

## ---- A1: frm_bootstrap() ------------------------------------------
set.seed(7)
sims <- simulate(fit, nsim = 40, re_formula = NA)
tops <- vapply(sims, function(v) max(as.integer(v)), 1L)
say("A1: simulate() nsim=40 seed 7: replicates whose max < 4: ",
    sum(tops < 4L), " of 40")

bs <- tryCatch(frm_bootstrap(fit, nsim = 40, seed = 7),
               error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(bs)) {
  say("A1: frm_bootstrap -> ", bs)
} else {
  say("A1: frm_bootstrap t dim = ", paste(dim(bs$t), collapse = "x"),
      "; colnames = ", paste(colnames(bs$t), collapse = ", "))
  say("A1: rows all-NA = ", sum(apply(bs$t, 1, function(r) all(is.na(r)))))
  say("A1: not converged = ", sum(!bs$converged))
  ## which replicates lacked the top category, and what came back
  bad <- which(tops < 4L)
  for (b in bad) {
    say("A1: replicate ", b, " max(y)=", tops[b], " t = ",
        paste(signif(bs$t[b, ], 6), collapse = " "))
  }
}

## ---- A2: refit() directly on a replicate that lost the category ---
bad <- which(tops < 4L)
if (length(bad)) {
  b <- bad[1L]
  yb <- as.integer(sims[[b]])
  rf <- tryCatch(suppressWarnings(refit(fit, yb)),
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(rf)) {
    say("A2: refit -> ", rf)
  } else {
    say("A2: refit thresholds = ", length(rf$estimates$tau_raw),
        " logLik = ", format(as.numeric(logLik(rf)), digits = 15))
  }
  ## the pinned comparison: a fresh frm() with thres(K) on that data
  db <- data.frame(x = dd$x, y = yb)
  fp <- tryCatch(suppressWarnings(
    frm(bf(y | thres(3) ~ x), family = cumulative(), data = db)),
    error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(fp)) {
    say("A2: thres(3) pinned fit -> ", fp)
  } else {
    say("A2: thres(3) pinned logLik = ",
        format(as.numeric(logLik(fp)), digits = 15),
        " thresholds = ", length(fp$estimates$tau_raw))
  }
  ## and the default (recounting) fresh fit
  fd <- tryCatch(suppressWarnings(
    frm(bf(y ~ x), family = cumulative(), data = db)),
    error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(fd)) {
    say("A2: default fresh fit -> ", fd)
  } else {
    say("A2: default fresh logLik = ",
        format(as.numeric(logLik(fd)), digits = 15),
        " thresholds = ", length(fd$estimates$tau_raw))
  }
}

## ---- A3: influence() leave-one-out --------------------------------
## the top category is reached by few rows; deleting one of them can
## remove it entirely
n4 <- sum(dd$y == 4L)
say("A3: rows with y == 4: ", n4)
d1 <- dd
if (n4 > 1L) {
  ## force the single-row case so the deletion is decisive
  keep <- which(dd$y == 4L)[1L]
  d1$y[setdiff(which(dd$y == 4L), keep)] <- 3L
}
say("A3: table(y) after forcing one top row = ",
    paste(table(d1$y), collapse = "/"))
f1 <- frm(bf(y ~ x), family = cumulative(), data = d1)
say("A3: fit thresholds = ", length(f1$estimates$tau_raw))
inf <- tryCatch(suppressWarnings(influence(f1, force = TRUE)),
                error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(inf)) {
  say("A3: influence -> ", inf)
} else {
  say("A3: influence fixed dim = ", paste(dim(inf$fixed), collapse = "x"))
  say("A3: colnames = ", paste(colnames(inf$fixed), collapse = ", "))
  ## the row for the deleted top-category observation
  itop <- which(d1$y == 4L)
  say("A3: deleted-top row = ", itop, " values = ",
      paste(signif(inf$fixed[itop, ], 6), collapse = " "))
  say("A3: any NA in that row = ", anyNA(inf$fixed[itop, ]))
  cd <- cooks.distance(inf)
  say("A3: cooks.distance at that row = ", signif(cd[itop], 6),
      "; max elsewhere = ", signif(max(cd[-itop], na.rm = TRUE), 6))
}

## ------------------------------------------------------------------
## B: a prior draw on an UNORDERED threshold vector
## ------------------------------------------------------------------
ddb <- data.frame(x = rnorm(40), y = 0)
for (fam in c("cratio", "acat", "sratio", "cumulative")) {
  f <- get(fam, envir = asNamespace("frmtmb"))()
  r <- tryCatch({
    s <- frm_simulate(bf(y | thres(3) ~ x), family = f, data = ddb,
                      prior = set_prior("normal(0, 2)",
                                        class = "Intercept") +
                        set_prior("normal(0, 1)", class = "b"),
                      nsim = 3, seed = 5)
    list(ok = TRUE, pars = attr(s, "pars"), y = s)
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  if (!r$ok) {
    say("B: ", fam, " -> ERROR: ", r$msg)
  } else {
    say("B: ", fam, " -> pars columns: ",
        paste(names(r$pars), collapse = ", "))
    print(r$pars)
    say("B: ", fam, " simulated category tables:")
    print(vapply(r$y, function(v) table(factor(v, levels = 1:4)),
                 integer(4)))
  }
}
