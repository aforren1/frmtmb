# Reviewer, punch round 1: B2, cmc = FALSE on mm() against brms standata.
# Z_k Z_k' (n x n) is compared, which does not depend on level order.
.libPaths(c("C:/Users/adf44/source/r/wt-formula2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages(expr)),
           error = function(e) structure(conditionMessage(e), class = "ERR"))
}
set.seed(31)
n <- 200
d <- data.frame(f = factor(sample(c("a", "b", "c", "d"), n, TRUE)),
                x = rnorm(n),
                h = factor(sample(sprintf("L%02d", 1:10), n, TRUE)),
                h2 = factor(sample(sprintf("L%02d", 1:10), n, TRUE)),
                w1 = runif(n, 0.2, 1), w2 = runif(n, 0.2, 1))
d$y <- rnorm(n, as.numeric(d$f) + d$x)
B <- brms::bf
cmp <- function(label, lhs, mmargs, cmc) {
  rhs <- sprintf("y ~ x + (%s | mm(%s))", lhs, mmargs)
  ff <- bf(as.formula(rhs), cmc = cmc)
  fb <- B(as.formula(rhs), cmc = cmc)
  fr <- q(frm(ff, data = d, dry_run = "frame"))
  sd <- q(brms::standata(fb, d))
  if (inherits(fr, "ERR") || inherits(sd, "ERR")) {
    cat(label, "cmc", cmc, "frm:", if (inherits(fr, "ERR")) fr else "ok",
        "| brms:", if (inherits(sd, "ERR")) sd else "ok", "\n")
    return(invisible())
  }
  cn <- unlist(lapply(fr$re_blocks, function(b)
    lapply(b$components, `[[`, "cnms")))
  K <- length(cn)
  Z <- as.matrix(fr$linpreds[["y.mu"]]$Z)
  M <- length(grep("^W_1_[0-9]+$", names(sd)))
  L <- sd$N_1
  kb <- length(grep("^Z_1_[0-9]+_1$", names(sd)))
  ok <- K == kb
  if (ok) for (k in seq_len(K)) {
    Zf <- Z[, seq(k, ncol(Z), by = K), drop = FALSE]
    Zb <- matrix(0, n, L)
    for (m in seq_len(M)) {
      J <- sd[[paste0("J_1_", m)]]
      v <- sd[[paste0("W_1_", m)]] * sd[[paste0("Z_1_", k, "_", m)]]
      Zb[cbind(seq_len(n), J)] <- Zb[cbind(seq_len(n), J)] + v
    }
    ok <- ok && isTRUE(all.equal(unname(tcrossprod(Zf)),
                                 unname(tcrossprod(Zb))))
  }
  cat(sprintf("%-26s cmc=%-5s frm cnms [%s] brms coefs %d  Z Z' equal %s\n",
              label, cmc, paste(cn, collapse = ","), kb, ok))
}
for (cmc in c(TRUE, FALSE)) {
  cmp("0 + f, equal weights", "0 + f", "h, h2", cmc)
  cmp("0 + f, weights", "0 + f", "h, h2, weights = cbind(w1, w2)", cmc)
  cmp("0 + f:x", "0 + f:x", "h, h2", cmc)
  cmp("0 + f + x, weights", "0 + f + x", "h, h2, weights = cbind(w1, w2)",
      cmc)
  cmp("0 + x", "0 + x", "h, h2", cmc)
  cmp("1 + f", "1 + f", "h, h2", cmc)
}

cat("== a fit and prediction\n")
fit <- frm(bf(y ~ x + (0 + f | mm(h, h2, weights = cbind(w1, w2))),
              cmc = FALSE), data = d)
print(variables(fit))
nd <- d[1:8, ]
cat("frm_linpred newdata == in-sample:",
    isTRUE(all.equal(unname(frm_linpred(fit, newdata = nd)),
                     unname(frm_linpred(fit)[1:8]))), "\n")
nd_new <- nd; nd_new$h <- factor(c("NEW", as.character(nd$h[-1])))
p1 <- q(fitted(fit, newdata = nd_new, allow_new_levels = TRUE))
cat("unseen member level, allow_new_levels:",
    if (inherits(p1, "ERR")) p1 else "ok", "\n")
p2 <- q(fitted(fit, newdata = nd_new))
cat("unseen member level, default:", if (inherits(p2, "ERR"))
  substr(p2, 1, 120) else "ok (no refusal)", "\n")
nd_f <- nd; nd_f$f <- factor(c("e", as.character(nd$f[-1])))
p3 <- q(fitted(fit, newdata = nd_f))
cat("unseen factor level e in the term:", if (inherits(p3, "ERR"))
  substr(p3, 1, 120) else "ok (no refusal)", "\n")
nd_a <- nd[nd$f == "a", ][1:2, ]; nd_a$f <- factor(nd_a$f)
p4 <- q(fitted(fit, newdata = nd_a))
cat("newdata with only level a:", if (inherits(p4, "ERR")) p4 else "ok",
    "\n")
p5 <- q(fitted(fit, newdata = nd,
               re_formula = ~ (0 + f | mm(h, h2, weights = cbind(w1, w2)))))
cat("re_formula naming the mm term:",
    if (inherits(p5, "ERR")) substr(p5, 1, 160) else
      isTRUE(all.equal(p5, fitted(fit, newdata = nd))), "\n")
s1 <- q(simulate(fit, nsim = 1, seed = 1, newdata = nd))
cat("simulate(newdata):", if (inherits(s1, "ERR")) s1 else "ok", "\n")
r <- q(ranef(fit)); cat("ranef cols:", colnames(r[[1]]), "\n")
cat("DONE\n")
