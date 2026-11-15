# Re-check of the refits fix round 2. ARM=before reads rellib-r21 and runs
# only part 4b (influence of non-refit fits, saved for a bitwise compare).
source("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-lib.R")
out <- "C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews"
q <- function(e) suppressWarnings(e)
tab <- function(i) list(fixed = unname(as.matrix(i$fixed)),
                        theta = unname(as.matrix(i$theta)))
mx <- function(a, b) max(abs(a - b), na.rm = TRUE)

cat("\n==== 4b influence of non-refit fits (both arms) ====\n")
set.seed(41)
n <- 60
db <- data.frame(x = rnorm(n), g = factor(rep(1:10, each = 6)),
                 w = runif(n, 0.5, 2), t = rep(1:6, 10))
db$y <- 1 + 0.5 * db$x + rnorm(10)[db$g] + rnorm(n)
db$s <- rbinom(n, 1, plogis(db$x))
db$keep <- db$x > -1.2
db$k <- rbinom(n, 8, plogis(0.2 + db$x))
db$ym <- db$y; db$ym[c(4, 30)] <- NA
db <- db[sample(n), ]                       # rows not in time order
nr <- list(
  re_groups = function() influence(frm(y ~ x + (1 | g), data = db), groups = "g"),
  bern12 = function() influence(frm(y12 ~ x, data = transform(db, y12 = s + 1),
                                    family = bernoulli()), force = TRUE),
  weights = function() influence(frm(y | weights(w) ~ x, data = db), force = TRUE),
  subset = function() influence(frm(y | subset(keep) ~ x, data = db), force = TRUE),
  trials = function() influence(frm(k | trials(8) ~ x, data = db,
                                    family = binomial()), force = TRUE),
  ar = function() influence(frm(y ~ x + ar(time = t, gr = g), data = db),
                            groups = "g"),
  logy = function() influence(frm(log(y + 10) ~ x, data = db), force = TRUE),
  mi = function() influence(frm(bf(ym | mi() ~ x), data = db), force = TRUE)
)
res <- lapply(names(nr), function(k) {
  r <- tryCatch(tab(q(nr[[k]]())), error = function(e) conditionMessage(e))
  cat(k, ": ", if (is.character(r)) paste("ERROR", substr(r, 1, 160)) else
        paste(dim(r$fixed), collapse = " x "), "\n", sep = "")
  r
})
names(res) <- names(nr)
saveRDS(res, file.path(out, paste0("refits-short-r2-infl-", arm, ".rds")))
if (identical(arm, "before")) quit(save = "no")

cat("\n==== 1 M1 / X1 / U1 U2: see refits-short-r2-extra.log ====\n")

cat("\n==== 2 ignore_response does not accept changed non-response data ====\n")
G <- 20
mkd <- function() {
  set.seed(1)
  dd <- data.frame(g = factor(rep(seq_len(G), each = 6)),
                   cl = factor(rep(seq_len(G / 2), each = 12)),
                   x = rnorm(G * 6), z = rnorm(G * 6))
  dd$y <- 1 + 0.5 * dd$x + rnorm(G * 6) + rnorm(G / 2)[dd$cl]
  dd
}
se_of <- function(v) round(sqrt(diag(v)), 7)
dd <- mkd()
fit <- frm(y ~ x, data = dd)
set.seed(2)
rf <- refit(fit, rnorm(nrow(dd)))
go("N1a refit, data unchanged: ~ cl = explicit",
   list(se_of(vcov_cluster(rf, ~ cl)), se_of(vcov_cluster(rf, dd$cl))))
dd$x[5] <- dd$x[5] + 1
go("N1b refit, dd$x changed: vcov_cluster(rf, ~ cl)", se_of(vcov_cluster(rf, ~ cl)))
go("N1c fit, dd$x changed: vcov_cluster(fit, ~ cl)", se_of(vcov_cluster(fit, ~ cl)))
go("N1d fit, dd$x changed: update(fit, . ~ . + z)", fixef(update(fit, . ~ . + z)))
dd <- mkd(); dd <- dd[c(2, 1, 3:nrow(dd)), ]; rownames(dd) <- NULL
go("N1e refit, two rows swapped, names reset: ~ cl", se_of(vcov_cluster(rf, ~ cl)))
dd <- mkd(); dd$y[5] <- dd$y[5] + 1
go("N1f fit, only dd$y changed: vcov_cluster(fit, ~ cl) (rows and cl intact)",
   list(se_of(vcov_cluster(fit, ~ cl)), se_of(vcov_cluster(fit, mkd()$cl))))
go("N1g fit, only dd$y changed: update(fit, . ~ . + z) still refused",
   fixef(update(fit, . ~ . + z)))
dd <- mkd(); dd <- dd[-7, ]
go("N1h fit, a row deleted: ~ cl", se_of(vcov_cluster(fit, ~ cl)))
# transformed response: the raw variable is still compared
dd <- mkd()
fl <- frm(log(y + 10) ~ x, data = dd)
dd$y[5] <- dd$y[5] + 1
go("N1i log(y + 10) fit, dd$y changed: ~ cl (raw y still compared?)",
   se_of(vcov_cluster(fl, ~ cl)))
# matrix response: a changed predictor is refused, a changed trials too
set.seed(5)
nm <- 200
xm <- rnorm(nm)
p <- cbind(a = 1, b = exp(0.3 + 0.4 * xm), c = exp(-0.2 * xm))
p <- p / rowSums(p)
Y <- t(vapply(seq_len(nm), function(i) rmultinom(1, 12, p[i, ]), numeric(3)))
colnames(Y) <- c("a", "b", "c")
dm <- data.frame(x = xm, N = 12L, cl = factor(rep(1:20, each = 10)))
dm$Y <- Y
fm <- frm(Y | trials(N) ~ x, data = dm, family = multinomial())
set.seed(6)
Y2 <- Y[sample(nm), ]
rm2 <- refit(fm, Y2)
go("N1j matrix refit, data unchanged: ~ cl = explicit",
   list(se_of(vcov_cluster(rm2, ~ cl)), se_of(vcov_cluster(rm2, dm$cl))))
dm0 <- dm
dm$x[3] <- dm$x[3] + 1
go("N1k matrix refit, dm$x changed: ~ cl", se_of(vcov_cluster(rm2, ~ cl)))
dm <- dm0; dm$N[3] <- 13L
go("N1l matrix refit, dm$N changed: ~ cl", se_of(vcov_cluster(rm2, ~ cl)))
dm <- dm0

cat("\n==== 3 N2 does not refuse legitimate fallbacks ====\n")
dd <- mkd()
ftop <- frm(y ~ x, data = dd)
go("N2a top-level fit, data present: ~ cl = explicit",
   list(se_of(vcov_cluster(ftop, ~ cl)), se_of(vcov_cluster(ftop, dd$cl))))
cl_ws <- factor(rep(1:5, each = 24))
go("N2b top-level fit, data present but lack cl_ws: ~ cl_ws = explicit",
   list(se_of(vcov_cluster(ftop, ~ cl_ws)), se_of(vcov_cluster(ftop, cl_ws))))
yv <- dd$y; xv <- dd$x; clv <- dd$cl
fnull <- go("N2c fit with data = NULL, variables global",
            class(frm(yv ~ xv, data = NULL)))
if (!is.character(fnull)) {
  fnull <- frm(yv ~ xv, data = NULL)
  go("N2c data = NULL, variables global: ~ clv = explicit",
     list(se_of(vcov_cluster(fnull, ~ clv)), se_of(vcov_cluster(fnull, clv))))
  go("N2d data = NULL: getCall()[['data']]", getCall(fnull)[["data"]])
}
fenv <- go("N2c' fit with data = globalenv()",
           class(frm(yv ~ xv, data = globalenv())))
if (!is.character(fenv)) {
  fenv <- frm(yv ~ xv, data = globalenv())
  go("N2c' data = globalenv(): ~ clv = explicit",
     list(se_of(vcov_cluster(fenv, ~ clv)), se_of(vcov_cluster(fenv, clv))))
}
fval <- do.call(frm, list(y ~ x, data = dd))
go("N2e data passed by value (do.call): ~ cl = explicit",
   list(se_of(vcov_cluster(fval, ~ cl)), se_of(vcov_cluster(fval, dd$cl))))
mkfit <- function(d) frm(y ~ x, data = d)
fa <- mkfit(dd)
go("N2f fit in a function on its argument (gone): ~ cl",
   se_of(vcov_cluster(fa, ~ cl)))
fd2 <- frm(y ~ x, data = dd, data2 = list(A = diag(2)))
go("N2g fit with data2 (call$data not data2): ~ cl = explicit",
   list(se_of(vcov_cluster(fd2, ~ cl)), se_of(vcov_cluster(fd2, dd$cl))))
go("N2h ~ cl inside a function, data global",
   (function(f) se_of(vcov_cluster(f, ~ cl)))(ftop))

cat("\n==== 4a influence of refits: equal to a fresh fit on the new response ====\n")
set.seed(8)
nn <- 30
du <- data.frame(x = rnorm(nn), w = runif(nn, 0.5, 2), N = 10L,
                 g = factor(rep(1:6, each = 5)), t = rep(1:5, 6))
du$y <- rnorm(nn)
du$k <- rbinom(nn, 10, 0.4)
du$s <- rbinom(nn, 1, 0.5)
du <- du[sample(nn), ]
ynew <- rnorm(nn, 2)
knew <- rbinom(nn, 10, 0.7)
snew <- rbinom(nn, 1, 0.3)
cmp <- function(lab, mkf, newresp, col, newcol, ...) {
  f <- mkf(du)
  if (is.function(newresp)) newresp <- newresp(f)
  r <- tryCatch(q(refit(f, newresp)), error = function(e) conditionMessage(e))
  if (is.character(r)) return(go(lab, paste("REFIT ERROR", r)))
  d2 <- du; d2[[col]] <- newcol
  a <- tryCatch(tab(q(influence(r, ...))), error = function(e) conditionMessage(e))
  b <- tab(q(influence(mkf(d2), ...)))
  o <- tab(q(influence(f, ...)))
  if (is.character(a)) return(go(lab, paste("ERROR", a)))
  go(lab, c(vs_fresh = mx(a$fixed, b$fixed), vs_original = mx(a$fixed, o$fixed),
            theta_vs_fresh = if (length(a$theta)) mx(a$theta, b$theta) else NA))
}
cmp("R1 gaussian", function(d) frm(y ~ x, data = d), ynew, "y", ynew, force = TRUE)
cmp("R2 weights(w)", function(d) frm(y | weights(w) ~ x, data = d), ynew, "y",
    ynew, force = TRUE)
cmp("R3 trials(N) binomial", function(d) frm(k | trials(N) ~ x, data = d,
                                             family = binomial()),
    knew, "k", knew, force = TRUE)
cmp("R4 bernoulli 1/2, refit given values", function(d)
  frm(y12 ~ x, data = transform(d, y12 = s + 1), family = bernoulli()),
  snew + 1, "s", snew, force = TRUE)
cmp("R5 (1 | g) groups = g", function(d) frm(y ~ x + (1 | g), data = d), ynew,
    "y", ynew, groups = "g")
cmp("R6 ar(time = t, gr = g), rows not in time order",
    function(d) frm(y ~ x + ar(time = t, gr = g), data = d), ynew, "y", ynew,
    groups = "g")
cmp("R7 log(y + 10) response", function(d) frm(log(y + 10) ~ x, data = d),
    log(ynew + 10), "y", ynew, force = TRUE)
# matrix response
set.seed(9)
Ym <- t(rmultinom(nn, 10, c(0.2, 0.3, 0.5))); colnames(Ym) <- c("a", "b", "c")
Ym2 <- t(rmultinom(nn, 10, c(0.5, 0.3, 0.2))); colnames(Ym2) <- colnames(Ym)
du$Ym <- Ym
cmp("R8 multinomial matrix response", function(d)
  frm(Ym | trials(N) ~ x, data = d, family = multinomial()), Ym2, "Ym", Ym2,
  force = TRUE)
# subset(): the refit's rows are the kept rows
du$keep <- du$x > -1
fs0 <- frm(y | subset(keep) ~ x, data = du)
cat("R9 subset fit: frame y length", length(fs0$frame$y[[1]]), "of", nn,
    "rows,", sum(du$keep), "kept\n")
ys <- if (length(fs0$frame$y[[1]]) == nn) ynew else ynew[du$keep]
cmp("R9 subset(keep)", function(d) frm(y | subset(keep) ~ x, data = d),
    ys, "y", if (length(ys) == nn) ynew else replace(du$y, du$keep, ys),
    force = TRUE)
# NA rows dropped at the fit
du$xna <- du$x; du$xna[c(2, 9)] <- NA
cmp("R10 na.action dropped 2 rows", function(d) frm(y ~ xna, data = d),
    ynew[-c(2, 9)], "y", replace(du$y, -c(2, 9), ynew[-c(2, 9)]), force = TRUE)
# multivariate fit: refit refuses it, so influence() keeps the old path
fmv <- frm(bf(y ~ x) + bf(k ~ x), data = du)
go("R11 refit of an mvbf fit", refit(fmv, ynew))
go("R11 influence(mvbf) runs", dim(q(influence(fmv, force = TRUE))$fixed))
# influence(refit, data = ) given explicitly: the given data win
f1 <- frm(y ~ x, data = du); r1 <- refit(f1, ynew)
go("R12 influence(refit, data = du) = influence(f1) (data given: theirs)",
   mx(tab(q(influence(r1, data = du, force = TRUE)))$fixed,
      tab(q(influence(f1, force = TRUE)))$fixed))
